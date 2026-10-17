# Helpers shared across analysis scripts.

suppressPackageStartupMessages(library(survey))

# Replicate weight matrix aligned to the rows of `dt` (n x 1000).
rep_weights <- function(dt, bsw) {
  idx <- match(dt$ADM_RNO, bsw$ADM_RNO)
  stopifnot(!anyNA(idx))
  as.matrix(bsw[idx, paste0("BSW", seq_len(n_boot)), with = FALSE])
}

# Bootstrap replicate design for model fitting with the survey package.
# Replicate weights are applied as supplied (combined with the sampling
# weight). Variance is computed around the full-sample estimate,
# v = (1/B) * sum((theta_b - theta)^2), B = 1000.
make_design <- function(dt, bsw) {
  svrepdesign(
    data = as.data.frame(dt), weights = ~WTS_M, repweights = rep_weights(dt, bsw),
    type = "bootstrap", combined.weights = TRUE,
    scale = 1 / n_boot, rscales = rep(1, n_boot), mse = TRUE
  )
}

# Weighted ratios sum(w * num) / sum(w * den) for each column of the 0/1
# matrices `num` and `den`, with the same bootstrap variance as make_design().
# Equivalent to survey::svymean on a domain, computed in one matrix product.
boot_ratio <- function(num, den, w, rw) {
  est <- colSums(num * w) / colSums(den * w)
  reps <- crossprod(rw, num) / crossprod(rw, den)
  se <- sqrt(colMeans(sweep(reps, 2, est)^2))
  data.table(p = est, se = se, n_num = colSums(num), n_den = colSums(den))
}

# Release rules from the CCHS 2017-2018 PUMF Quality Assurance Guidelines:
# unweighted count < 30 -> not released (S); CV 15.1-35.0% -> E; CV > 35.0% -> F.
release_flag <- function(n_unweighted, cv) {
  fcase(
    n_unweighted < 30, "S",
    cv > 0.35,         "F",
    cv > 0.15,         "E",
    default = ""
  )
}

# Weighted prevalence of a 0/1 variable within domains defined by `by`
# (character vector of factor columns, or NULL for the whole sample), with a
# 95% CI on the logit scale (delta method) and release flags.
prevalence <- function(dt, rw, var, by = NULL) {
  y <- dt[[var]]
  dom <- if (is.null(by)) data.table(all = factor("All")) else dt[, ..by]
  key <- interaction(dom, drop = TRUE, sep = "|")
  lev <- levels(key)
  den <- sapply(lev, function(l) as.numeric(key == l & !is.na(y) & !is.na(key)))
  num <- den * (y %in% 1L)
  res <- boot_ratio(num, den, dt$WTS_M, rw)
  res <- cbind(as.data.table(tstrsplit(lev, "|", fixed = TRUE)), res)
  setnames(res, seq_along(names(dom)), names(dom))
  res[, `:=`(
    variable = var,
    cv = se / p,
    lcl = plogis(qlogis(p) - 1.96 * se / (p * (1 - p))),
    ucl = plogis(qlogis(p) + 1.96 * se / (p * (1 - p)))
  )]
  res[, flag := release_flag(pmin(n_num, n_den), cv)]
  res[]
}

# Weighted distribution of a factor within domains (missing kept as a level).
distribution <- function(dt, rw, var, by) {
  x <- addNA(dt[[var]], ifany = TRUE)
  levels(x)[is.na(levels(x))] <- "Missing"
  key <- interaction(dt[, ..by], drop = TRUE, sep = "|")
  combos <- CJ(dom = levels(key), lev = levels(x), sorted = FALSE)
  den <- sapply(combos$dom, function(l) as.numeric(key == l & !is.na(key)))
  num <- den * sapply(combos$lev, function(l) as.numeric(x == l))
  res <- boot_ratio(num, den, dt$WTS_M, rw)
  res <- cbind(combos, res)
  res[, (by) := tstrsplit(dom, "|", fixed = TRUE)]
  res[, `:=`(variable = var, level = lev, cv = se / p)]
  res[, flag := release_flag(n_num, cv)]
  res[, c("dom", "lev") := NULL]
  res[]
}

fmt_pct <- function(p, digits = 1) formatC(100 * p, format = "f", digits = digits)

fmt_prev <- function(p, lcl, ucl, flag) {
  fifelse(flag %in% c("S", "F"), flag,
          paste0(fmt_pct(p), " (", fmt_pct(lcl), "-", fmt_pct(ucl), ")",
                 fifelse(flag == "E", " E", "")))
}

# Logistic regression ----------------------------------------------------------

# Weighted logistic regression by Newton-Raphson. Weights are rescaled to
# mean 1 by the caller (as survey::svyglm does), which keeps the iterations
# numerically stable with population-scale weights.
irls_logit <- function(X, y, w, beta = rep(0, ncol(X)), tol = 1e-8, maxit = 50L) {
  for (i in seq_len(maxit)) {
    mu <- plogis(drop(X %*% beta))
    step <- drop(solve(crossprod(X, X * (w * mu * (1 - mu))), crossprod(X, w * (y - mu))))
    beta <- beta + step
    if (max(abs(step)) < tol) return(list(beta = beta, iter = i, converged = TRUE))
  }
  list(beta = beta, iter = maxit, converged = FALSE)
}

# Survey-weighted logistic regression with bootstrap replicate variance.
# All predictors are categorical, so records sharing a covariate pattern and
# outcome are collapsed and their full-sample and replicate weights summed;
# the weighted likelihood, and so every estimate, is unchanged.
# Returns the full-sample coefficients, the B x p matrix of replicate
# coefficients and the bootstrap covariance (1/B) * sum((b_r - b)(b_r - b)').
fit_logit_boot <- function(dt, rw, outcome, rhs) {
  vars <- unique(c(all.vars(rhs), outcome))
  g <- dt[, .(grp = .GRP), by = vars]
  grp <- g[dt, on = vars, x.grp]
  stopifnot(!anyNA(grp))
  w_full <- rowsum(dt$WTS_M, grp, reorder = TRUE)[, 1]
  w_full <- w_full / (sum(dt$WTS_M) / nrow(dt))
  w_rep  <- rowsum(rw, grp, reorder = TRUE)
  w_rep  <- sweep(w_rep, 2, colSums(rw) / nrow(dt), "/")
  setorder(g, grp)
  X <- model.matrix(rhs, g)
  y <- g[[outcome]]
  full <- irls_logit(X, y, w_full)
  if (!full$converged) stop("Full-sample fit did not converge: ", outcome)
  reps <- t(vapply(seq_len(ncol(w_rep)), function(r) {
    f <- irls_logit(X, y, w_rep[, r], full$beta)
    if (!f$converged) rep(NA_real_, ncol(X)) else f$beta
  }, numeric(ncol(X))))
  colnames(reps) <- colnames(X)
  # Replicates whose fit did not converge are dropped from the variance and
  # counted in n_fail.
  ok <- complete.cases(reps)
  dev <- sweep(reps[ok, , drop = FALSE], 2, full$beta)
  list(beta = setNames(full$beta, colnames(X)), reps = reps,
       vcov = crossprod(dev, dev) / sum(ok), n_fail = sum(!ok),
       n = nrow(dt), n_patterns = nrow(g))
}

# Multiple imputation pooling ------------------------------------------------

# Rubin's rules for a scalar: q and u are the per-imputation estimates and
# within-imputation variances. With one data set (complete-case analysis) the
# result is the complete-data estimate with a normal reference distribution.
rubin <- function(q, u) {
  m <- length(q)
  if (anyNA(q) || anyNA(u)) {
    return(data.table(est = if (anyNA(q)) NA_real_ else mean(q), se = NA_real_, df = NA_real_,
                      lcl = NA_real_, ucl = NA_real_, p = NA_real_, fmi = NA_real_))
  }
  qbar <- mean(q); ubar <- mean(u)
  b <- if (m > 1) var(q) else 0
  tv <- ubar + (1 + 1 / m) * b
  r <- (1 + 1 / m) * b / ubar
  df <- if (m > 1 && b > 0) (m - 1) * (1 + 1 / r)^2 else Inf
  se <- sqrt(tv)
  data.table(est = qbar, se = se, df = df,
             lcl = qbar - qt(0.975, df) * se, ucl = qbar + qt(0.975, df) * se,
             p = 2 * pt(abs(qbar / se), df, lower.tail = FALSE),
             fmi = if (m > 1) (r + 2 / (df + 3)) / (r + 1) else 0)
}

# Rubin's rules for linear combinations L %*% beta. `fits` is a list (one
# element per imputation) of objects with $beta and $vcov. L is a k x p
# matrix with column names matching the coefficients.
pool_contrasts <- function(fits, L) {
  L <- L[, names(fits[[1]]$beta), drop = FALSE]
  rbindlist(lapply(seq_len(nrow(L)), function(i) {
    l <- L[i, ]
    q <- vapply(fits, function(f) sum(l * f$beta), numeric(1))
    u <- vapply(fits, function(f) drop(t(l) %*% f$vcov %*% l), numeric(1))
    rubin(q, u)[, term := rownames(L)[i]]
  }))[, c("term", "est", "se", "df", "lcl", "ucl", "p", "fmi")]
}

# Pool a nonlinear function of the coefficients. `fun` maps a matrix of
# coefficients (rows = full sample or replicates, named columns) to a named
# list or data.frame of measures. Within each imputation the variance is the
# bootstrap variance of the measure across replicates around its full-sample
# value; imputations are combined with Rubin's rules.
pool_function <- function(fits, fun) {
  per_imp <- lapply(fits, function(f) {
    full <- as.data.table(fun(matrix(f$beta, nrow = 1, dimnames = list(NULL, names(f$beta)))))
    reps <- as.data.table(fun(f$reps[complete.cases(f$reps), , drop = FALSE]))
    v <- vapply(names(full), function(nm) mean((reps[[nm]] - full[[nm]])^2), numeric(1))
    list(q = unlist(full), u = v)
  })
  nms <- names(per_imp[[1]]$q)
  rbindlist(lapply(nms, function(nm) {
    rubin(vapply(per_imp, function(x) x$q[[nm]], numeric(1)),
          vapply(per_imp, function(x) x$u[[nm]], numeric(1)))[, measure := nm]
  }))[, c("measure", "est", "se", "df", "lcl", "ucl", "p", "fmi")]
}

# Additive interaction between BMI category `b` and sex from a model with a
# bmi_cat x sex product term (references: normal weight, male). Following
# Knol et al. (2011) the reference is the lowest-risk group, normal-weight
# women, so the second exposure is male sex:
#   OR10 = overweight/obese women, OR01 = normal-weight men,
#   OR11 = overweight/obese men, each against normal-weight women.
# RERI = OR11 - OR10 - OR01 + 1; AP = RERI / OR11;
# SI = (OR11 - 1) / ((OR10 - 1) + (OR01 - 1)), pooled on the log scale.
additive_or <- function(b) {
  bt <- paste0("bmi_cat", b); it <- paste0("bmi_cat", b, ":sexFemale")
  function(B) {
    or10 <- exp(B[, bt] + B[, it])
    or01 <- exp(-B[, "sexFemale"])
    or11 <- exp(B[, bt] - B[, "sexFemale"])
    reri <- or11 - or10 - or01 + 1
    # SI is undefined when its numerator or denominator is not positive; if
    # that happens in any replicate its variance is left undefined.
    si <- (or11 - 1) / ((or10 - 1) + (or01 - 1))
    list(log_or10 = log(or10), log_or01 = log(or01), log_or11 = log(or11),
         reri = reri, ap = reri / or11,
         log_si = ifelse(si > 0, log(pmax(si, .Machine$double.xmin)), NA_real_))
  }
}

# Model-standardized prevalence for each sex and BMI category, averaging
# predicted probabilities over the weighted covariate distribution of the
# whole cohort (marginal standardization). Returns, for BMI category `b`
# against normal weight, the standardized prevalence differences in women and
# men and their difference (women minus men), i.e. the interaction contrast
# on the probability scale. Replicate estimates use the replicate weights for
# both the coefficients and the standardizing distribution. `rw` holds the
# replicate weights for the records identified by `rw_ids` (ADM_RNO); each
# data set may be any subset of them.
std_prevalence <- function(fits, data_list, rw, rw_ids, rhs, b, by_vars) {
  per_imp <- lapply(seq_along(fits), function(k) {
    f <- fits[[k]]; dt <- data_list[[k]]
    idx <- match(dt$ADM_RNO, rw_ids)
    stopifnot(!anyNA(idx))
    pat <- dt[, .(grp = .GRP), by = by_vars]
    grp <- pat[dt, on = by_vars, x.grp]
    w_full <- rowsum(dt$WTS_M, grp, reorder = TRUE)[, 1]
    w_rep <- rowsum(rw[idx, , drop = FALSE], grp, reorder = TRUE)
    setorder(pat, grp)
    scen <- CJ(sex = c("Female", "Male"), bmi = c("Normal", b), sorted = FALSE)
    pr <- lapply(seq_len(nrow(scen)), function(i) {
      nd <- copy(pat)
      nd[, sex := factor(scen$sex[i], levels(dt$sex))]
      nd[, bmi_cat := factor(scen$bmi[i], levels(dt$bmi_cat))]
      X <- model.matrix(rhs, nd)[, names(f$beta), drop = FALSE]
      ok <- complete.cases(f$reps)
      full <- sum(w_full * plogis(drop(X %*% f$beta))) / sum(w_full)
      reps <- colSums(w_rep[, ok, drop = FALSE] * plogis(X %*% t(f$reps[ok, , drop = FALSE]))) /
        colSums(w_rep[, ok, drop = FALSE])
      list(full = full, reps = reps)
    })
    names(pr) <- paste(scen$sex, scen$bmi, sep = "_")
    get_m <- function(x) {
      p_fn <- x[[paste0("Female_Normal")]]; p_fb <- x[[paste0("Female_", b)]]
      p_mn <- x[[paste0("Male_Normal")]];   p_mb <- x[[paste0("Male_", b)]]
      list(p_female_normal = p_fn, p_female_bmi = p_fb, p_male_normal = p_mn, p_male_bmi = p_mb,
           rd_female = p_fb - p_fn, rd_male = p_mb - p_mn, rd_diff = (p_fb - p_fn) - (p_mb - p_mn))
    }
    full <- get_m(lapply(pr, `[[`, "full"))
    reps <- get_m(lapply(pr, `[[`, "reps"))
    list(q = unlist(full), u = vapply(names(full), function(nm) mean((reps[[nm]] - full[[nm]])^2), numeric(1)))
  })
  nms <- names(per_imp[[1]]$q)
  rbindlist(lapply(nms, function(nm) {
    rubin(vapply(per_imp, function(x) x$q[[nm]], numeric(1)),
          vapply(per_imp, function(x) x$u[[nm]], numeric(1)))[, measure := nm]
  }))[, c("measure", "est", "se", "df", "lcl", "ucl", "p", "fmi")]
}

# Multivariate Wald test pooled across imputations (D1; Li, Raghunathan and
# Rubin 1991) for H0: beta[terms] = 0.
pool_wald <- function(fits, terms) {
  m <- length(fits); k <- length(terms)
  q <- sapply(fits, function(f) f$beta[terms])
  if (k == 1L) q <- matrix(q, nrow = 1)
  qbar <- rowMeans(q)
  ubar <- Reduce(`+`, lapply(fits, function(f) f$vcov[terms, terms, drop = FALSE])) / m
  b <- if (m == 1) matrix(0, k, k) else if (k == 1L) matrix(var(q[1, ])) else cov(t(q))
  r <- (1 + 1 / m) * sum(diag(b %*% solve(ubar))) / k
  d1 <- drop(t(qbar) %*% solve(ubar) %*% qbar) / (k * (1 + r))
  t_ <- k * (m - 1)
  v <- if (m == 1) Inf else if (t_ > 4) 4 + (t_ - 4) * (1 + (1 - 2 / t_) / r)^2 else t_ * (1 + 1 / k) * (1 + 1 / r)^2 / 2
  data.table(statistic = d1, df1 = k, df2 = v, p = pf(d1, k, v, lower.tail = FALSE))
}

fmt_or <- function(est, lcl, ucl, digits = 2) {
  f <- function(x) formatC(exp(x), format = "f", digits = digits)
  paste0(f(est), " (", f(lcl), "-", f(ucl), ")")
}

fmt_p <- function(p) fifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3))

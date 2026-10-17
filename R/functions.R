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
  n_fail <- sum(!complete.cases(reps))
  dev <- sweep(reps, 2, full$beta)
  list(beta = setNames(full$beta, colnames(X)), reps = reps,
       vcov = crossprod(dev, dev) / nrow(reps), n_fail = n_fail,
       n = nrow(dt), n_patterns = nrow(g))
}

# Multiple imputation pooling ------------------------------------------------

# Rubin's rules for linear combinations L %*% beta. `fits` is a list (one
# element per imputation) of objects with $beta and $vcov. L is a k x p
# matrix with column names matching the coefficients.
pool_contrasts <- function(fits, L) {
  m <- length(fits)
  L <- L[, names(fits[[1]]$beta), drop = FALSE]
  q <- vapply(fits, function(f) drop(L %*% f$beta), numeric(nrow(L)))
  u <- vapply(fits, function(f) diag(L %*% f$vcov %*% t(L)), numeric(nrow(L)))
  if (is.null(dim(q))) { q <- matrix(q, nrow = 1); u <- matrix(u, nrow = 1) }
  qbar <- rowMeans(q); ubar <- rowMeans(u); b <- apply(q, 1, var)
  tv <- ubar + (1 + 1 / m) * b
  r <- (1 + 1 / m) * b / ubar
  df <- (m - 1) * (1 + 1 / r)^2
  se <- sqrt(tv)
  data.table(term = rownames(L), est = qbar, se = se, df = df,
             lcl = qbar - qt(0.975, df) * se, ucl = qbar + qt(0.975, df) * se,
             p = 2 * pt(abs(qbar / se), df, lower.tail = FALSE),
             fmi = (r + 2 / (df + 3)) / (r + 1))
}

# Multivariate Wald test pooled across imputations (D1; Li, Raghunathan and
# Rubin 1991) for H0: beta[terms] = 0.
pool_wald <- function(fits, terms) {
  m <- length(fits); k <- length(terms)
  q <- sapply(fits, function(f) f$beta[terms])
  if (k == 1L) q <- matrix(q, nrow = 1)
  qbar <- rowMeans(q)
  ubar <- Reduce(`+`, lapply(fits, function(f) f$vcov[terms, terms, drop = FALSE])) / m
  b <- if (k == 1L) matrix(var(q[1, ])) else cov(t(q))
  r <- (1 + 1 / m) * sum(diag(b %*% solve(ubar))) / k
  d1 <- drop(t(qbar) %*% solve(ubar) %*% qbar) / (k * (1 + r))
  t_ <- k * (m - 1)
  v <- if (t_ > 4) 4 + (t_ - 4) * (1 + (1 - 2 / t_) / r)^2 else t_ * (1 + 1 / k) * (1 + 1 / r)^2 / 2
  data.table(statistic = d1, df1 = k, df2 = v, p = pf(d1, k, v, lower.tail = FALSE))
}

fmt_or <- function(est, lcl, ucl, digits = 2) {
  f <- function(x) formatC(exp(x), format = "f", digits = digits)
  paste0(f(est), " (", f(lcl), "-", f(ucl), ")")
}

fmt_p <- function(p) fifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3))

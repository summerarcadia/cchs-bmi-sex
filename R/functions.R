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

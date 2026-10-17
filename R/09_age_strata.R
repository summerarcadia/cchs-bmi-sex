# Sex x BMI interaction within age groups (18-44, 45-64, 65+).
#
# Outcomes: multimorbidity (primary) and diabetes (the outcome with a
# multiplicative interaction in the whole cohort). Underweight adults are
# excluded: there are 0-1 multimorbidity cases among underweight adults aged
# 18-44, which makes the category inestimable. Exclusion uses each
# imputation's BMI. Models adjust for age (finer groups within the stratum),
# household income, household education, smoking, alcohol and physical
# activity.
#
# Set REUSE_FITS=1 to rebuild the tables from previously saved fits.
#
# Outputs (not tracked): data/derived/age_strata_fits.rds
# Tracked: output/tables/age_strata_or.csv, output/tables/age_strata_additive.csv,
#          output/tables/age_strata_tests.csv

source("R/00_setup.R")
source("R/functions.R")
suppressPackageStartupMessages(library(parallel))

n_core <- max(1L, detectCores() - 2L)

imputed <- readRDS(file.path(derived_dir, "imputed.rds"))
ids <- imputed[[1]]$ADM_RNO
prep <- function(d) {
  d <- d[bmi_cat != "Underweight"]
  d[, `:=`(bmi_cat = factor(as.character(bmi_cat), c("Normal", "Overweight", "Obese")),
           sex = relevel(sex, "Male"))]
  d[]
}
bsw <- readRDS(file.path(derived_dir, "bsw.rds"))
rw  <- rep_weights(imputed[[1]], bsw)
rm(bsw); invisible(gc())

outcomes <- c(cmm = "Multimorbidity", dm = "Diabetes")
strata <- c("18-44", "45-64", "65+")
covars <- c("age_adj", "income", "educ", "smoking", "alcohol", "phys_act")
rhs <- as.formula(paste("~ bmi_cat * sex +", paste(covars, collapse = " + ")))
specs <- CJ(outcome = names(outcomes), stratum = strata, sorted = FALSE)

fits_path <- file.path(derived_dir, "age_strata_fits.rds")
if (Sys.getenv("REUSE_FITS") == "1" && file.exists(fits_path)) {
  fits <- readRDS(fits_path)
} else {
  t0 <- Sys.time()
  fits <- mclapply(seq_along(imputed), function(k) {
    d_all <- prep(imputed[[k]])
    lapply(strata, function(s) {
      d <- d_all[age_grp == s]
      d[, age_adj := droplevels(age_adj)]
      rwk <- rw[match(d$ADM_RNO, ids), , drop = FALSE]
      lapply(names(outcomes), function(o) fit_logit_boot(d, rwk, o, rhs))
    })
  }, mc.cores = n_core, mc.preschedule = FALSE)
  cat("Model fitting time:", format(Sys.time() - t0), "\n")
  err <- vapply(fits, inherits, logical(1), "try-error")
  if (any(err)) stop("Model fitting failed: ", paste(unique(unlist(fits[err])), collapse = "; "))

  # fits[[outcome]][[stratum]][[imputation]]
  fits <- setNames(lapply(seq_along(outcomes), function(j) {
    setNames(lapply(seq_along(strata), function(i) lapply(fits, function(f) f[[i]][[j]])), strata)
  }), names(outcomes))
  n_fail <- sum(unlist(lapply(fits, function(o) lapply(o, function(s) vapply(s, `[[`, numeric(1), "n_fail")))))
  cat("Replicate fits that did not converge (dropped from variance):", n_fail, "\n")
  saveRDS(fits, fits_path)
}

# Unweighted case counts (observed data, underweight excluded).
co <- readRDS(file.path(derived_dir, "cohort.rds"))
cases <- rbindlist(lapply(names(outcomes), function(o) {
  co[!is.na(bmi_cat) & bmi_cat != "Underweight" & !is.na(get(o)),
     .(n_cases = sum(get(o) == 1L)), by = .(age_grp, sex, bmi_cat)][, outcome := o]
}))
min_cases <- cases[, .(min_cases = min(n_cases)), by = .(outcome, stratum = as.character(age_grp))]

# Odds ratios within each sex and ratio of ORs.
bmi_nonref <- c("Overweight", "Obese")
or_res <- rbindlist(lapply(names(outcomes), function(o) {
  rbindlist(lapply(strata, function(s) {
    bn <- names(fits[[o]][[s]][[1]]$beta)
    L <- do.call(rbind, c(
      setNames(lapply(bmi_nonref, function(b) setNames(as.numeric(bn == paste0("bmi_cat", b)), bn)), paste0("men|", bmi_nonref)),
      setNames(lapply(bmi_nonref, function(b) setNames(as.numeric(bn %in% c(paste0("bmi_cat", b), paste0("bmi_cat", b, ":sexFemale"))), bn)), paste0("women|", bmi_nonref)),
      setNames(lapply(bmi_nonref, function(b) setNames(as.numeric(bn == paste0("bmi_cat", b, ":sexFemale")), bn)), paste0("ratio|", bmi_nonref))
    ))
    pool_contrasts(fits[[o]][[s]], L)[, `:=`(outcome = o, stratum = s)]
  }))
}))
or_res[, c("contrast", "bmi") := tstrsplit(term, "|", fixed = TRUE)]
or_res[min_cases, on = .(outcome, stratum), min_cases := i.min_cases]

tests <- rbindlist(lapply(names(outcomes), function(o) {
  rbindlist(lapply(strata, function(s) {
    pool_wald(fits[[o]][[s]], paste0("bmi_cat", bmi_nonref, ":sexFemale"))[, `:=`(outcome = o, stratum = s)]
  }))
}))
tests[min_cases, on = .(outcome, stratum), min_cases := i.min_cases]

# Additive interaction and standardized prevalence differences.
add_res <- rbindlist(lapply(names(outcomes), function(o) {
  rbindlist(lapply(strata, function(s) {
    rbindlist(lapply(bmi_nonref, function(b) {
      a <- pool_function(fits[[o]][[s]], additive_or(b))
      data_s <- lapply(imputed, function(d) { x <- prep(d)[age_grp == s]; x[, age_adj := droplevels(age_adj)]; x })
      sp <- std_prevalence(fits[[o]][[s]], data_s, rw, ids, rhs, b, covars)
      rbind(a, sp)[, `:=`(outcome = o, stratum = s, bmi = b)]
    }))
  }))
}))
add_res[min_cases, on = .(outcome, stratum), min_cases := i.min_cases]
add_res[, `:=`(
  value = fifelse(startsWith(measure, "log_"), exp(est), est),
  value_lcl = fifelse(startsWith(measure, "log_"), exp(lcl), lcl),
  value_ucl = fifelse(startsWith(measure, "log_"), exp(ucl), ucl)
)]
add_res[, measure := sub("^log_", "", measure)]

fwrite(or_res[, .(outcome, stratum, contrast, bmi, or = round(exp(est), 3), lcl = round(exp(lcl), 3),
                  ucl = round(exp(ucl), 3), p, df = round(df, 1), fmi = round(fmi, 3), min_cases, sparse = min_cases < 30)],
       file.path(tab_dir, "age_strata_or.csv"))
fwrite(tests[, .(outcome, stratum, statistic = round(statistic, 3), df1, df2 = round(df2, 1), p, min_cases, sparse = min_cases < 30)],
       file.path(tab_dir, "age_strata_tests.csv"))
fwrite(add_res[, .(outcome, stratum, bmi, measure, value = round(value, 4), lcl = round(value_lcl, 4),
                   ucl = round(value_ucl, 4), p, df = round(df, 1), fmi = round(fmi, 3), min_cases, sparse = min_cases < 30)],
       file.path(tab_dir, "age_strata_additive.csv"))

print(tests)
print(add_res[measure %in% c("reri", "ap", "rd_female", "rd_male", "rd_diff"),
              .(outcome, stratum, bmi, measure, est = round(value, 3), lcl = round(value_lcl, 3),
                ucl = round(value_ucl, 3), p = fmt_p(p), sparse = min_cases < 30)])

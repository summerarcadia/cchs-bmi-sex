# Sensitivity analyses for the primary outcome (multimorbidity).
#
#   complete_case     adults with no missing data (no imputation)
#   excl_underweight  imputed data, underweight adults excluded
#   self_report       BMI category from self-reported BMI (no bias adjustment)
#   obesity_classes   adjusted BMI with obesity split into class I (30-34.9)
#                     and class II-III (>= 35)
#
# Each repeats the interaction model of R/06_models.R and the additive
# interaction measures of R/08_additive.R.
#
# Set REUSE_FITS=1 to rebuild the tables from previously saved fits.
#
# Outputs (not tracked): data/derived/sensitivity_fits.rds
# Tracked: output/tables/sensitivity.csv, output/tables/sensitivity_tests.csv

source("R/00_setup.R")
source("R/functions.R")
suppressPackageStartupMessages(library(parallel))

n_core <- max(1L, detectCores() - 2L)

imputed <- readRDS(file.path(derived_dir, "imputed.rds"))
ids <- imputed[[1]]$ADM_RNO
bsw <- readRDS(file.path(derived_dir, "bsw.rds"))
rw  <- rep_weights(imputed[[1]], bsw)
rm(bsw); invisible(gc())
co <- readRDS(file.path(derived_dir, "cohort.rds"))

covars <- c("age_adj", "income", "educ", "smoking", "alcohol", "phys_act")
rhs <- as.formula(paste("~ bmi_cat * sex +", paste(covars, collapse = " + ")))

set_levels <- function(d, levs) {
  d[, `:=`(bmi_cat = factor(as.character(bmi_cat), levs), sex = relevel(sex, "Male"))]
  stopifnot(!anyNA(d$bmi_cat))
  d[]
}
std_levels <- c("Normal", "Underweight", "Overweight", "Obese")

variants <- list(
  complete_case = list(
    data = list(set_levels(co[complete == TRUE], std_levels)),
    levels = c("Underweight", "Overweight", "Obese")
  ),
  excl_underweight = list(
    data = lapply(imputed, function(d) set_levels(copy(d)[bmi_cat != "Underweight"], c("Normal", "Overweight", "Obese"))),
    levels = c("Overweight", "Obese")
  ),
  self_report = list(
    data = lapply(imputed, function(d) { d <- copy(d); d[, bmi_cat := bmi_cat_sr]; set_levels(d, std_levels) }),
    levels = c("Underweight", "Overweight", "Obese")
  ),
  obesity_classes = list(
    data = lapply(imputed, function(d) {
      d <- copy(d)
      d[, bmi_cat := cut(bmi_adj, c(-Inf, 18.5, 25, 30, 35, Inf),
                         c("Underweight", "Normal", "Overweight", "ObeseI", "ObeseII"), right = FALSE)]
      set_levels(d, c("Normal", "Underweight", "Overweight", "ObeseI", "ObeseII"))
    }),
    levels = c("Underweight", "Overweight", "ObeseI", "ObeseII")
  )
)
rm(imputed); invisible(gc())

# Fit -------------------------------------------------------------------------

fits_path <- file.path(derived_dir, "sensitivity_fits.rds")
if (Sys.getenv("REUSE_FITS") == "1" && file.exists(fits_path)) {
  fits <- readRDS(fits_path)
} else {
  t0 <- Sys.time()
  fits <- lapply(names(variants), function(v) {
    dl <- variants[[v]]$data
    res <- mclapply(seq_along(dl), function(k) {
      d <- dl[[k]]
      fit_logit_boot(d, rw[match(d$ADM_RNO, ids), , drop = FALSE], "cmm", rhs)
    }, mc.cores = n_core, mc.preschedule = FALSE)
    err <- vapply(res, inherits, logical(1), "try-error")
    if (any(err)) stop("Fitting failed in ", v, ": ", paste(unique(unlist(res[err])), collapse = "; "))
    res
  })
  names(fits) <- names(variants)
  cat("Model fitting time:", format(Sys.time() - t0), "\n")
  print(vapply(fits, function(f) sum(vapply(f, `[[`, numeric(1), "n_fail")), numeric(1)))
  saveRDS(fits, fits_path)
}

# Summarise -------------------------------------------------------------------

summarise_variant <- function(v) {
  f <- fits[[v]]; levs <- variants[[v]]$levels; dl <- variants[[v]]$data
  bn <- names(f[[1]]$beta)
  unit <- function(terms) setNames(as.numeric(bn %in% terms), bn)
  L <- do.call(rbind, c(
    setNames(lapply(levs, function(b) unit(paste0("bmi_cat", b))), paste0("men|", levs)),
    setNames(lapply(levs, function(b) unit(c(paste0("bmi_cat", b), paste0("bmi_cat", b, ":sexFemale")))), paste0("women|", levs)),
    setNames(lapply(levs, function(b) unit(paste0("bmi_cat", b, ":sexFemale"))), paste0("ratio|", levs))
  ))
  ors <- pool_contrasts(f, L)
  ors[, c("measure", "bmi") := tstrsplit(term, "|", fixed = TRUE)]
  ors[, measure := paste0("or_", measure)]
  ors[, `:=`(value = exp(est), value_lcl = exp(lcl), value_ucl = exp(ucl))]
  add_levs <- setdiff(levs, "Underweight")
  add <- rbindlist(lapply(add_levs, function(b) {
    a <- pool_function(f, additive_or(b))[measure %in% c("reri", "ap", "log_si")]
    sp <- std_prevalence(f, dl, rw, ids, rhs, b, covars)[measure %in% c("rd_female", "rd_male", "rd_diff")]
    rbind(a, sp)[, bmi := b]
  }))
  add[, `:=`(
    value = fifelse(measure == "log_si", exp(est), est),
    value_lcl = fifelse(measure == "log_si", exp(lcl), lcl),
    value_ucl = fifelse(measure == "log_si", exp(ucl), ucl)
  )]
  add[measure == "log_si", measure := "si"]
  test <- pool_wald(f, paste0("bmi_cat", levs, ":sexFemale"))
  list(
    est = rbind(ors[, .(bmi, measure, value, value_lcl, value_ucl, p, fmi)],
                add[, .(bmi, measure, value, value_lcl, value_ucl, p, fmi)])[, `:=`(variant = v, n = f[[1]]$n, m = length(f))],
    test = test[, `:=`(variant = v, n = f[[1]]$n, m = length(f))]
  )
}
out <- lapply(names(variants), summarise_variant)
est <- rbindlist(lapply(out, `[[`, "est"))
tests <- rbindlist(lapply(out, `[[`, "test"))

fwrite(est[, .(variant, n, m, bmi, measure, value = round(value, 4), lcl = round(value_lcl, 4),
               ucl = round(value_ucl, 4), p, fmi = round(fmi, 3))],
       file.path(tab_dir, "sensitivity.csv"))
fwrite(tests[, .(variant, n, m, statistic = round(statistic, 3), df1, df2 = round(df2, 1), p)],
       file.path(tab_dir, "sensitivity_tests.csv"))

print(tests)
print(est[bmi != "Underweight" & measure %in% c("or_men", "or_women", "or_ratio", "reri", "rd_diff"),
          .(variant, bmi, measure, est = round(value, 3), lcl = round(value_lcl, 3), ucl = round(value_ucl, 3), p = fmt_p(p))])

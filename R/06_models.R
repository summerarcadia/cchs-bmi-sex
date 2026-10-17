# Survey-weighted logistic regression on the multiply imputed data.
#
# Outcomes: cardiometabolic multimorbidity (primary); hypertension, diabetes,
# heart disease (secondary). For each outcome two models are fitted:
#   main:        BMI category + sex + covariates
#   interaction: BMI category x sex + covariates
# Covariates: age group, household income, household education, smoking,
# alcohol, physical activity. Reference levels: normal weight, male.
#
# Variance: 1,000 bootstrap replicate weights within each imputation, pooled
# across 20 imputations with Rubin's rules.
#
# Output (not tracked): data/derived/model_fits.rds
# Tables are built from the saved fits in R/07_model_tables.R.

source("R/00_setup.R")
source("R/functions.R")
suppressPackageStartupMessages(library(parallel))

n_core <- max(1L, detectCores() - 2L)

imputed <- readRDS(file.path(derived_dir, "imputed.rds"))
imputed <- imputed[seq_len(as.integer(Sys.getenv("MODELS_M", length(imputed))))]
m <- length(imputed)
for (k in seq_len(m)) {
  imputed[[k]][, `:=`(bmi_cat = relevel(bmi_cat, "Normal"), sex = relevel(sex, "Male"))]
}
stopifnot(all(vapply(imputed, function(x) identical(x$ADM_RNO, imputed[[1]]$ADM_RNO), logical(1))))

bsw <- readRDS(file.path(derived_dir, "bsw.rds"))
rw  <- rep_weights(imputed[[1]], bsw)
rm(bsw); invisible(gc())

outcomes <- c(cmm = "Multimorbidity", htn = "Hypertension", dm = "Diabetes", hd = "Heart disease")
covars <- "age_adj + income + educ + smoking + alcohol + phys_act"
rhs <- list(
  main        = as.formula(paste("~ bmi_cat + sex +", covars)),
  interaction = as.formula(paste("~ bmi_cat * sex +", covars))
)
specs <- CJ(outcome = names(outcomes), model = names(rhs), sorted = FALSE)

# Fit -------------------------------------------------------------------------

t0 <- Sys.time()
fits <- mclapply(seq_len(m), function(k) {
  lapply(seq_len(nrow(specs)), function(j) {
    f <- fit_logit_boot(imputed[[k]], rw, specs$outcome[j], rhs[[specs$model[j]]])
    if (specs$model[j] == "main") f$reps <- NULL
    f
  })
}, mc.cores = n_core, mc.preschedule = FALSE)
cat("Model fitting time:", format(Sys.time() - t0), "\n")
err <- vapply(fits, inherits, logical(1), "try-error")
if (any(err)) stop("Model fitting failed: ", paste(unique(unlist(fits[err])), collapse = "; "))

# Reorganise as fits[[outcome]][[model]][[imputation]].
fits <- setNames(lapply(names(outcomes), function(o) {
  setNames(lapply(names(rhs), function(md) {
    j <- which(specs$outcome == o & specs$model == md)
    lapply(fits, `[[`, j)
  }), names(rhs))
}), names(outcomes))

n_fail <- sum(unlist(lapply(fits, function(o) lapply(o, function(md) vapply(md, `[[`, numeric(1), "n_fail")))))
cat("Replicate fits that did not converge:", n_fail, "\n")
stopifnot(n_fail == 0)
saveRDS(fits, file.path(derived_dir, "model_fits.rds"))

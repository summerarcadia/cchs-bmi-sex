# Propensity score weighting for BMI category, within sex.
#
# Within each sex, a multinomial logistic propensity model for BMI category
# (normal, overweight, obese; underweight excluded) on age group, household
# income, household education, smoking, alcohol and physical activity, with
# the survey weight as sampling weight (WeightIt, method "glm", estimand ATE,
# stabilized). Outcomes: multimorbidity and diabetes, weighted by survey
# weight x propensity weight, in two forms:
#   marginal       BMI category x sex only. Each sex's BMI contrast refers to
#                  that sex's own covariate distribution, and the comparison
#                  between sexes is not covariate-adjusted.
#   doubly_robust  BMI category x sex plus the adjustment covariates, with
#                  prevalence standardized to the whole cohort, as in the main
#                  analysis.
#
# Variance: the replicate weights are multiplied by the propensity weights,
# which are held fixed (not re-estimated in each replicate), so the CIs do
# not reflect uncertainty in the propensity model.
#
# Set REUSE_FITS=1 to rebuild the outputs from previously saved fits.
#
# Outputs (not tracked): data/derived/propensity_fits.rds
# Tracked: output/tables/propensity.csv, output/tables/propensity_balance.csv,
#          output/figures/balance_men.png, output/figures/balance_women.png

source("R/00_setup.R")
source("R/functions.R")
suppressPackageStartupMessages({
  library(WeightIt)
  library(cobalt)
  library(ggplot2)
  library(parallel)
})

# Each worker holds a propensity-weighted copy of the replicate weights.
n_core <- min(4L, max(1L, detectCores() - 2L))

imputed <- readRDS(file.path(derived_dir, "imputed.rds"))
ids <- imputed[[1]]$ADM_RNO
bsw <- readRDS(file.path(derived_dir, "bsw.rds"))
rw  <- rep_weights(imputed[[1]], bsw)
rm(bsw); invisible(gc())

covars <- c("age_adj", "income", "educ", "smoking", "alcohol", "phys_act")
ps_formula <- as.formula(paste("bmi_cat ~", paste(covars, collapse = " + ")))
outcomes <- c(cmm = "Multimorbidity", dm = "Diabetes")
rhs <- list(
  marginal = ~ bmi_cat * sex,
  doubly_robust = as.formula(paste("~ bmi_cat * sex +", paste(covars, collapse = " + ")))
)
approaches <- names(rhs)
bmi_nonref <- c("Overweight", "Obese")

cov_labels <- c(
  age_adj = "Age", income = "Income", educ = "Education", smoking = "Smoking",
  alcohol = "Alcohol", phys_act = "Physical activity"
)

prep <- function(d) {
  d <- copy(d)[bmi_cat != "Underweight"]
  d[, `:=`(bmi_cat = factor(as.character(bmi_cat), c("Normal", "Overweight", "Obese")),
           sex = relevel(sex, "Male"))]
  d[]
}

fits_path <- file.path(derived_dir, "propensity_fits.rds")
if (Sys.getenv("REUSE_FITS") == "1" && file.exists(fits_path)) {
  per_imp <- readRDS(fits_path)
} else {
  per_imp <- mclapply(seq_along(imputed), function(k) {
    d <- prep(imputed[[k]])
    d[, psw := NA_real_]
    bal <- list()
    for (s in levels(d$sex)) {
      ds <- d[sex == s]
      W <- weightit(ps_formula, data = ds, method = "glm", estimand = "ATE",
                    s.weights = ds$WTS_M, stabilize = TRUE)
      d[sex == s, psw := W$weights]
      bt <- bal.tab(W, stats = "m", binary = "std", un = TRUE, thresholds = c(m = 0.1))
      ess <- as.data.table(bt$Observations, keep.rownames = "type")
      bal[[s]] <- list(
        max_smd = data.table(sex = s, imputation = k,
                             max_smd_unadjusted = max(abs(bt$Balance.Across.Pairs$Max.Diff.Un)),
                             max_smd_weighted = max(abs(bt$Balance.Across.Pairs$Max.Diff.Adj)),
                             ess_normal = ess[2, Normal], ess_overweight = ess[2, Overweight], ess_obese = ess[2, Obese],
                             psw_max = max(W$weights), psw_p99 = quantile(W$weights, 0.99)),
        W = if (k == 1L) W else NULL
      )
    }
    stopifnot(!anyNA(d$psw))
    dw <- copy(d)[, WTS_M := WTS_M * psw]
    rwk <- rw[match(d$ADM_RNO, ids), , drop = FALSE] * d$psw
    f <- lapply(approaches, function(a) {
      setNames(lapply(names(outcomes), function(o) fit_logit_boot(dw, rwk, o, rhs[[a]])), names(outcomes))
    })
    names(f) <- approaches
    rm(rwk); invisible(gc())
    list(fits = f, data = d, balance = rbindlist(lapply(bal, `[[`, "max_smd")),
         W = lapply(bal, `[[`, "W"))
  }, mc.cores = n_core, mc.preschedule = FALSE)
  err <- vapply(per_imp, inherits, logical(1), "try-error")
  if (any(err)) stop("Propensity analysis failed: ", paste(unique(unlist(per_imp[err])), collapse = "; "))
  saveRDS(per_imp, fits_path)
}

balance <- rbindlist(lapply(per_imp, `[[`, "balance"))
fwrite(balance, file.path(tab_dir, "propensity_balance.csv"))

# Love plots from the first imputation.
for (s in c("Male", "Female")) {
  W <- per_imp[[1]]$W[[s]]
  bt <- bal.tab(W, stats = "m", binary = "std", un = TRUE)
  vars <- rownames(bt$Balance.Across.Pairs)
  pretty <- vars
  for (v in names(cov_labels)) pretty <- sub(paste0("^", v, "_"), paste0(cov_labels[[v]], ": "), pretty)
  p <- love.plot(W, stats = "m", binary = "std", abs = TRUE, thresholds = c(m = 0.1),
                 var.order = "unadjusted", var.names = setNames(pretty, vars),
                 colors = c("grey55", "black"),
                 sample.names = c("Survey-weighted", "Survey x propensity weighted"),
                 title = paste0("Covariate balance across BMI categories, ", ifelse(s == "Male", "men", "women")),
                 subtitle = "Maximum absolute standardized mean difference across category pairs (imputation 1)") +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, paste0("balance_", ifelse(s == "Male", "men", "women"), ".png")),
         p, width = 6.5, height = 6, dpi = 300, bg = "white")
}

# Pooled estimates.
# Standardization uses the survey-weighted covariate distribution of the
# whole (non-underweight) cohort; for the marginal model the predictions do
# not depend on covariates.
data_list <- lapply(per_imp, `[[`, "data")
std_by <- list(marginal = "sex", doubly_robust = covars)

res <- rbindlist(lapply(approaches, function(a) {
  rbindlist(lapply(names(outcomes), function(o) {
    fo <- lapply(per_imp, function(x) x$fits[[a]][[o]])
    bn <- names(fo[[1]]$beta)
    unit <- function(terms) setNames(as.numeric(bn %in% terms), bn)
    L <- do.call(rbind, c(
      setNames(lapply(bmi_nonref, function(b) unit(paste0("bmi_cat", b))), paste0("or_men|", bmi_nonref)),
      setNames(lapply(bmi_nonref, function(b) unit(c(paste0("bmi_cat", b), paste0("bmi_cat", b, ":sexFemale")))), paste0("or_women|", bmi_nonref)),
      setNames(lapply(bmi_nonref, function(b) unit(paste0("bmi_cat", b, ":sexFemale"))), paste0("or_ratio|", bmi_nonref))
    ))
    ors <- pool_contrasts(fo, L)
    ors[, c("measure", "bmi") := tstrsplit(term, "|", fixed = TRUE)]
    ors[, `:=`(value = exp(est), value_lcl = exp(lcl), value_ucl = exp(ucl))]
    add <- rbindlist(lapply(bmi_nonref, function(b) {
      ad <- pool_function(fo, additive_or(b))[measure %in% c("reri", "ap", "log_si")]
      sp <- std_prevalence(fo, data_list, rw, ids, rhs[[a]], b, std_by[[a]])[measure %in% c("rd_female", "rd_male", "rd_diff")]
      rbind(ad, sp)[, bmi := b]
    }))
    add[, `:=`(value = fifelse(measure == "log_si", exp(est), est),
               value_lcl = fifelse(measure == "log_si", exp(lcl), lcl),
               value_ucl = fifelse(measure == "log_si", exp(ucl), ucl))]
    add[measure == "log_si", measure := "si"]
    test <- pool_wald(fo, paste0("bmi_cat", bmi_nonref, ":sexFemale"))
    rbind(ors[, .(bmi, measure, value, value_lcl, value_ucl, p, fmi)],
          add[, .(bmi, measure, value, value_lcl, value_ucl, p, fmi)],
          data.table(bmi = "All", measure = "global_interaction_p", value = test$p, value_lcl = NA, value_ucl = NA, p = test$p, fmi = NA))[, `:=`(outcome = o, approach = a)]
  }))
}))
cat("Replicate fits that did not converge:",
    sum(unlist(lapply(per_imp, function(x) lapply(x$fits, function(a) vapply(a, `[[`, numeric(1), "n_fail"))))), "\n")
fwrite(res[, .(approach, outcome, bmi, measure, value = round(value, 4), lcl = round(value_lcl, 4),
               ucl = round(value_ucl, 4), p, fmi = round(fmi, 3))],
       file.path(tab_dir, "propensity.csv"))

print(balance[, .(max_smd_unadj = max(max_smd_unadjusted), max_smd_wtd = max(max_smd_weighted),
                  min_ess_obese = min(ess_obese), max_psw = max(psw_max)), by = sex])
print(res[measure %in% c("or_men", "or_women", "or_ratio", "reri", "rd_diff", "global_interaction_p"),
          .(approach, outcome, bmi, measure, est = round(value, 3), lcl = round(value_lcl, 3), ucl = round(value_ucl, 3), p = fmt_p(p))])

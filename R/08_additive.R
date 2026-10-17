# Additive interaction between BMI category and sex, from the interaction
# models fitted in R/06_models.R.
#
# Odds ratio scale: RERI, attributable proportion (AP) and synergy index (SI)
# with normal-weight women (the lowest-risk group) as the reference, so the
# second exposure is male sex (Knol et al. 2011).
# Probability scale: model-standardized prevalence of each outcome by sex and
# BMI category, the BMI-associated prevalence difference within each sex, and
# the difference between women and men (interaction contrast).
#
# The synergy index is reported without a CI when its denominator or
# numerator is not positive in any replicate (unstable when an OR is near 1).
#
# Variance: within each imputation, the bootstrap variance of each measure
# across the 1,000 replicate coefficient vectors (and, for standardized
# prevalence, the replicate weights); pooled with Rubin's rules.
#
# Outputs: output/tables/additive_interaction.csv,
#          output/tables/standardized_prevalence.csv

source("R/00_setup.R")
source("R/functions.R")

fits <- readRDS(file.path(derived_dir, "model_fits.rds"))
imputed <- readRDS(file.path(derived_dir, "imputed.rds"))
for (k in seq_along(imputed)) {
  imputed[[k]][, `:=`(bmi_cat = relevel(bmi_cat, "Normal"), sex = relevel(sex, "Male"))]
}
bsw <- readRDS(file.path(derived_dir, "bsw.rds"))
rw  <- rep_weights(imputed[[1]], bsw)
rm(bsw); invisible(gc())

outcomes <- c(cmm = "Multimorbidity", htn = "Hypertension", dm = "Diabetes", hd = "Heart disease")
covars <- c("age_adj", "income", "educ", "smoking", "alcohol", "phys_act")
rhs <- as.formula(paste("~ bmi_cat * sex +", paste(covars, collapse = " + ")))
bmi_levels <- c("Overweight", "Obese")

add <- rbindlist(lapply(names(outcomes), function(o) {
  rbindlist(lapply(bmi_levels, function(b) {
    pool_function(fits[[o]]$interaction, additive_or(b))[, `:=`(outcome = o, bmi = b)]
  }))
}))
add[, `:=`(
  value = fifelse(startsWith(measure, "log_"), exp(est), est),
  value_lcl = fifelse(startsWith(measure, "log_"), exp(lcl), lcl),
  value_ucl = fifelse(startsWith(measure, "log_"), exp(ucl), ucl)
)]
add[, measure := sub("^log_", "", measure)]
fwrite(add[, .(outcome, bmi, measure, value = round(value, 4), lcl = round(value_lcl, 4),
               ucl = round(value_ucl, 4), p, se_estimation_scale = se, df = round(df, 1), fmi = round(fmi, 3))],
       file.path(tab_dir, "additive_interaction.csv"))

std <- rbindlist(lapply(names(outcomes), function(o) {
  rbindlist(lapply(bmi_levels, function(b) {
    std_prevalence(fits[[o]]$interaction, imputed, rw, imputed[[1]]$ADM_RNO, rhs, b, covars)[, `:=`(outcome = o, bmi = b)]
  }))
}))
fwrite(std[, .(outcome, bmi, measure, pct = round(100 * est, 2), lcl_pct = round(100 * lcl, 2),
               ucl_pct = round(100 * ucl, 2), p, se_pct = round(100 * se, 3), df = round(df, 1), fmi = round(fmi, 3))],
       file.path(tab_dir, "standardized_prevalence.csv"))

print(add[measure %in% c("reri", "ap", "si"),
          .(outcome, bmi, measure, est = round(value, 3), lcl = round(value_lcl, 3), ucl = round(value_ucl, 3), p = fmt_p(p))])
print(std[measure %in% c("rd_female", "rd_male", "rd_diff"),
          .(outcome, bmi, measure, pct = round(100 * est, 2), lcl = round(100 * lcl, 2), ucl = round(100 * ucl, 2), p = fmt_p(p))])

# Formatted tables for the additive interaction, age-stratified, sensitivity,
# propensity score and E-value results (R/08 to R/12).
#
# Outputs: output/tables/table5.html, output/tables/tableS1_age_strata.html,
#          output/tables/tableS2_sensitivity.html, output/tables/tableS3_evalues.html

source("R/00_setup.R")
source("R/functions.R")
suppressPackageStartupMessages(library(gt))

outcomes <- c(cmm = "Multimorbidity", htn = "Hypertension", dm = "Diabetes", hd = "Heart disease")
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
ci <- function(v, l, u, d = 2) fifelse(is.na(l), paste0(f2(v, d), " (not estimable)"),
                                       paste0(f2(v, d), " (", f2(l, d), " to ", f2(u, d), ")"))
pp <- function(v, l, u) ci(100 * v, 100 * l, 100 * u, 1)
dagger <- function(x, flag) paste0(x, fifelse(flag, " &dagger;", ""))

# Table 5: additive interaction, whole cohort --------------------------------

add <- fread(file.path(tab_dir, "additive_interaction.csv"))
std <- fread(file.path(tab_dir, "standardized_prevalence.csv"))
std[, `:=`(value = pct / 100, lcl = lcl_pct / 100, ucl = ucl_pct / 100)]

t5 <- merge(
  dcast(add[measure %in% c("reri", "ap", "si")], outcome + bmi ~ measure,
        value.var = c("value", "lcl", "ucl")),
  dcast(std[measure %in% c("rd_female", "rd_male", "rd_diff")], outcome + bmi ~ measure,
        value.var = c("value", "lcl", "ucl")),
  by = c("outcome", "bmi")
)
t5 <- t5[, .(
  outcome = outcomes[outcome], bmi,
  rd_female = pp(value_rd_female, lcl_rd_female, ucl_rd_female),
  rd_male = pp(value_rd_male, lcl_rd_male, ucl_rd_male),
  rd_diff = pp(value_rd_diff, lcl_rd_diff, ucl_rd_diff),
  reri = ci(value_reri, lcl_reri, ucl_reri),
  ap = ci(value_ap, lcl_ap, ucl_ap),
  si = ci(value_si, lcl_si, ucl_si)
)]
t5 <- t5[order(match(outcome, outcomes), match(bmi, c("Overweight", "Obese")))]

t5_gt <- gt(t5, groupname_col = "outcome", rowname_col = "bmi") |>
  tab_header(title = "Table 5. Additive interaction between BMI category and sex",
             subtitle = "Each BMI category against normal weight. Estimate (95% CI)") |>
  tab_spanner("Standardized prevalence difference, percentage points", columns = c(rd_female, rd_male, rd_diff)) |>
  tab_spanner("Odds ratio scale (reference: normal-weight women)", columns = c(reri, ap, si)) |>
  cols_label(rd_female = "Women", rd_male = "Men", rd_diff = "Women minus men",
             reri = "RERI", ap = "AP", si = "SI") |>
  tab_source_note(paste(
    "Standardized prevalence: predicted prevalence from the sex x BMI interaction model averaged over the",
    "survey-weighted covariate distribution of the cohort. A negative difference (women minus men) means a",
    "larger absolute excess in men. RERI, AP and SI use normal-weight women, the lowest-risk group, as the",
    "reference, so the second exposure is male sex (Knol et al. 2011); RERI > 0, AP > 0 and SI > 1 indicate",
    "more than additive joint association of the BMI category and male sex. Odds ratios approximate risk ratios",
    "less well when an outcome is common (hypertension; multimorbidity at older ages). SI is not estimable when its",
    "denominator is not positive in some bootstrap replicates. Models adjusted for age group, household income,",
    "household education, smoking, alcohol and physical activity. Bootstrap variance (1,000 replicates) pooled",
    "over 20 imputations."
  ))
gtsave(t5_gt, file.path(tab_dir, "table5.html"))

# Table S1: age strata ----------------------------------------------------------

ao <- fread(file.path(tab_dir, "age_strata_or.csv"))
aa <- fread(file.path(tab_dir, "age_strata_additive.csv"))
at <- fread(file.path(tab_dir, "age_strata_tests.csv"))
s1 <- merge(
  dcast(ao[contrast %in% c("men", "women", "ratio")], outcome + stratum + bmi + sparse ~ contrast,
        value.var = c("or", "lcl", "ucl")),
  dcast(aa[measure %in% c("reri", "rd_diff")], outcome + stratum + bmi ~ measure,
        value.var = c("value", "lcl", "ucl")),
  by = c("outcome", "stratum", "bmi")
)
s1 <- merge(s1, at[, .(outcome, stratum, p_int = p)], by = c("outcome", "stratum"))
s1 <- s1[, .(
  group = paste0(outcomes[outcome], ", ", stratum, " (interaction p = ", fmt_p(p_int), ")"),
  bmi,
  men = dagger(paste0(f2(or_men), " (", f2(lcl_men), "-", f2(ucl_men), ")"), sparse),
  women = dagger(paste0(f2(or_women), " (", f2(lcl_women), "-", f2(ucl_women), ")"), sparse),
  ratio = dagger(paste0(f2(or_ratio), " (", f2(lcl_ratio), "-", f2(ucl_ratio), ")"), sparse),
  reri = dagger(ci(value_reri, lcl_reri, ucl_reri), sparse),
  rd_diff = dagger(pp(value_rd_diff, lcl_rd_diff, ucl_rd_diff), sparse),
  o = match(outcome, names(outcomes)), s = stratum
)]
setorder(s1, o, s, bmi)
s1[, bmi := factor(bmi, c("Overweight", "Obese"))]
setorder(s1, o, s, bmi)
s1[, c("o", "s") := NULL]

s1_gt <- gt(s1, groupname_col = "group", rowname_col = "bmi") |>
  tab_header(title = "Table S1. Sex x BMI interaction by age group",
             subtitle = "Underweight excluded. Estimate (95% CI)") |>
  tab_spanner("Odds ratio vs normal weight", columns = c(men, women, ratio)) |>
  cols_label(men = "Men", women = "Women", ratio = "Ratio (women / men)", reri = "RERI",
             rd_diff = "Prevalence difference, women minus men (pp)") |>
  fmt_markdown(columns = everything()) |>
  tab_source_note(html(paste(
    "Models fitted within each age group; adjustment and variance as in Table 5. Interaction p: pooled Wald test",
    "of the two product terms. &dagger; Fewer than 30 unweighted cases in at least one sex and BMI category of",
    "the stratum (multimorbidity at ages 18-44); not interpreted."
  )))
gtsave(s1_gt, file.path(tab_dir, "tableS1_age_strata.html"))

# Table S2: sensitivity analyses ---------------------------------------------

main_or <- fread(file.path(tab_dir, "models_interaction.csv"))[outcome == "cmm" & contrast %in% c("men", "women", "ratio") & level == "Obese"]
main_add <- add[outcome == "cmm" & bmi == "Obese" & measure == "reri"]
main_rd <- std[outcome == "cmm" & bmi == "Obese" & measure == "rd_diff"]
main_row <- data.table(
  variant = "Main analysis (adjusted BMI, multiple imputation)", bmi = "Obese",
  men = main_or[contrast == "men", fmt_or(log(or), log(lcl), log(ucl))],
  women = main_or[contrast == "women", fmt_or(log(or), log(lcl), log(ucl))],
  ratio = main_or[contrast == "ratio", fmt_or(log(or), log(lcl), log(ucl))],
  reri = ci(main_add$value, main_add$lcl, main_add$ucl),
  rd_diff = pp(main_rd$value, main_rd$lcl, main_rd$ucl)
)
sens <- fread(file.path(tab_dir, "sensitivity.csv"))
variant_labels <- c(
  complete_case = "Complete cases only",
  excl_underweight = "Underweight excluded",
  self_report = "Self-reported BMI (no bias adjustment)",
  obesity_classes = "Obesity class I (30-34.9) and II-III (>= 35)"
)
sens_rows <- sens[bmi %in% c("Obese", "ObeseI", "ObeseII")]
sw <- dcast(sens_rows[measure %in% c("or_men", "or_women", "or_ratio", "reri", "rd_diff")],
            variant + bmi ~ measure, value.var = c("value", "lcl", "ucl"))
sw <- sw[, .(
  variant = paste0(variant_labels[variant], " (n = ", format(sens[match(sw$variant, variant), n], big.mark = ","), ")"),
  bmi = fcase(bmi == "ObeseI", "Obese class I", bmi == "ObeseII", "Obese class II-III", default = bmi),
  men = fmt_or(log(value_or_men), log(lcl_or_men), log(ucl_or_men)),
  women = fmt_or(log(value_or_women), log(lcl_or_women), log(ucl_or_women)),
  ratio = fmt_or(log(value_or_ratio), log(lcl_or_ratio), log(ucl_or_ratio)),
  reri = ci(value_reri, lcl_reri, ucl_reri),
  rd_diff = pp(value_rd_diff, lcl_rd_diff, ucl_rd_diff),
  o = match(variant, names(variant_labels))
)]
setorder(sw, o, bmi); sw[, o := NULL]

ps <- fread(file.path(tab_dir, "propensity.csv"))[outcome == "cmm" & bmi == "Obese"]
ps_labels <- c(
  doubly_robust = "Propensity score weighted within sex, covariate-adjusted model (underweight excluded)",
  marginal = "Propensity score weighted within sex, marginal model (underweight excluded)"
)
psw <- dcast(ps[measure %in% c("or_men", "or_women", "or_ratio", "reri", "rd_diff")], approach + bmi ~ measure,
             value.var = c("value", "lcl", "ucl"))
ps_row <- psw[, .(
  variant = ps_labels[approach], bmi,
  men = fmt_or(log(value_or_men), log(lcl_or_men), log(ucl_or_men)),
  women = fmt_or(log(value_or_women), log(lcl_or_women), log(ucl_or_women)),
  ratio = fmt_or(log(value_or_ratio), log(lcl_or_ratio), log(ucl_or_ratio)),
  reri = ci(value_reri, lcl_reri, ucl_reri),
  rd_diff = pp(value_rd_diff, lcl_rd_diff, ucl_rd_diff),
  o = match(approach, names(ps_labels))
)]
setorder(ps_row, o); ps_row[, o := NULL]

s2 <- rbind(main_row, sw, ps_row)
s2_gt <- gt(s2, groupname_col = "variant", rowname_col = "bmi") |>
  tab_header(title = "Table S2. Sensitivity analyses: obesity, sex and multimorbidity",
             subtitle = "Obesity (or obesity class) against normal weight. Estimate (95% CI)") |>
  tab_spanner("Odds ratio vs normal weight", columns = c(men, women, ratio)) |>
  cols_label(men = "Men", women = "Women", ratio = "Ratio (women / men)", reri = "RERI",
             rd_diff = "Prevalence difference, women minus men (pp)") |>
  tab_source_note(paste(
    "RERI with normal-weight women as reference. All models covariate-adjusted except the marginal propensity score",
    "model, where covariates enter only through weights that balance them across BMI categories within each sex; its",
    "comparison between sexes is not covariate-adjusted and each sex's prevalence difference refers to that sex's own",
    "covariate distribution. Prevalence differences elsewhere are standardized to the whole cohort. Propensity weights",
    "are held fixed across bootstrap replicates, so those CIs do not include uncertainty from the propensity model."
  ))
gtsave(s2_gt, file.path(tab_dir, "tableS2_sensitivity.html"))

# Table S3: E-values -------------------------------------------------------------

ev <- fread(file.path(tab_dir, "evalues.csv"))
s3 <- ev[, .(outcome = outcomes[outcome], row = paste0(fifelse(sex == "men", "Men", "Women"), ", ", tolower(bmi)),
             or = fmt_or(log(or), log(lcl), log(ucl)), evalue = f2(evalue), evalue_ci = f2(evalue_ci),
             o = match(outcome, names(outcomes)))]
setorder(s3, o); s3[, o := NULL]
s3_gt <- gt(s3, groupname_col = "outcome", rowname_col = "row") |>
  tab_header(title = "Table S3. E-values for unmeasured confounding",
             subtitle = "Sex-specific odds ratios against normal weight") |>
  cols_label(or = "Odds ratio (95% CI)", evalue = "E-value, estimate", evalue_ci = "E-value, CI limit") |>
  tab_source_note(paste(
    "E-value: minimum strength of association, on the risk ratio scale, that an unmeasured confounder would need",
    "with both BMI category and the outcome, conditional on the measured covariates, to fully explain away the",
    "observed odds ratio (estimate) or move its CI to include 1 (CI limit). Odds ratios treated as risk ratios for",
    "outcomes with weighted prevalence below 15%; square-root conversion for hypertension (18.7%)."
  ))
gtsave(s3_gt, file.path(tab_dir, "tableS3_evalues.html"))
cat("Tables written\n")

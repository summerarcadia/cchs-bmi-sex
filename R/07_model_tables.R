# Odds ratio tables from the saved model fits (R/06_models.R).
#
# Multiple testing: Holm adjustment across the four global sex x BMI
# interaction tests (one per outcome); Benjamini-Hochberg q-values across the
# twelve category-specific ratio-of-odds-ratio tests.
#
# Outputs: output/tables/models_main.csv, output/tables/models_interaction.csv,
#          output/tables/interaction_tests.csv, output/tables/case_counts_sex_bmi.csv,
#          output/tables/table3.html, output/tables/table4.html

source("R/00_setup.R")
source("R/functions.R")
suppressPackageStartupMessages(library(gt))

fits <- readRDS(file.path(derived_dir, "model_fits.rds"))
outcomes <- c(cmm = "Multimorbidity", htn = "Hypertension", dm = "Diabetes", hd = "Heart disease")
stopifnot(identical(names(fits), names(outcomes)))
m <- length(fits$cmm$main)
cat("Imputations:", m, "\n")

# Unweighted case counts by sex and BMI (observed data) ----------------------

co <- readRDS(file.path(derived_dir, "cohort.rds"))
cases <- rbindlist(lapply(names(outcomes), function(o) {
  co[!is.na(bmi_cat) & !is.na(get(o)), .(n_cases = sum(get(o) == 1L), n = .N), by = .(sex, bmi_cat)][, outcome := o]
}))

# Main-effects models -------------------------------------------------------

coef_names <- names(fits$cmm$main[[1]]$beta)[-1]
L_main <- diag(length(coef_names) + 1)[-1, , drop = FALSE]
dimnames(L_main) <- list(coef_names, names(fits$cmm$main[[1]]$beta))

main_res <- rbindlist(lapply(names(outcomes), function(o) {
  pool_contrasts(fits[[o]]$main, L_main)[, outcome := o]
}))
main_res[, `:=`(or = exp(est), or_lcl = exp(lcl), or_ucl = exp(ucl))]
fwrite(main_res[, .(outcome, term, log_or = est, se, or = round(or, 3), lcl = round(or_lcl, 3),
                    ucl = round(or_ucl, 3), p, df = round(df, 1), fmi = round(fmi, 3))],
       file.path(tab_dir, "models_main.csv"))

# Interaction models ----------------------------------------------------------

bnames <- names(fits$cmm$interaction[[1]]$beta)
bmi_nonref <- c("Underweight", "Overweight", "Obese")
e <- function(x = numeric()) {
  v <- setNames(numeric(length(bnames)), bnames)
  stopifnot(all(names(x) %in% bnames))
  v[names(x)] <- v[names(x)] + x
  v
}
bt <- function(b) paste0("bmi_cat", b)
it <- function(b) paste0("bmi_cat", b, ":sexFemale")
L_int <- do.call(rbind, c(
  # BMI vs normal weight, within men and within women
  setNames(lapply(bmi_nonref, function(b) e(setNames(1, bt(b)))), paste0("men|", bmi_nonref)),
  setNames(lapply(bmi_nonref, function(b) e(setNames(c(1, 1), c(bt(b), it(b))))), paste0("women|", bmi_nonref)),
  # Ratio of odds ratios, women / men (multiplicative interaction)
  setNames(lapply(bmi_nonref, function(b) e(setNames(1, it(b)))), paste0("ratio|", bmi_nonref)),
  # Women vs men within each BMI category
  list("sex|Normal" = e(c(sexFemale = 1))),
  setNames(lapply(bmi_nonref, function(b) e(setNames(c(1, 1), c("sexFemale", it(b))))), paste0("sex|", bmi_nonref)),
  # Joint effects, single reference: normal-weight men
  list("joint_m|Male:Normal" = e()),
  setNames(lapply(bmi_nonref, function(b) e(setNames(1, bt(b)))), paste0("joint_m|Male:", bmi_nonref)),
  list("joint_m|Female:Normal" = e(c(sexFemale = 1))),
  setNames(lapply(bmi_nonref, function(b) e(setNames(c(1, 1, 1), c(bt(b), "sexFemale", it(b))))),
           paste0("joint_m|Female:", bmi_nonref))
))
stopifnot(nrow(L_int) == 21L, sum(abs(L_int)) == 3 + 6 + 3 + 1 + 6 + 0 + 3 + 1 + 9)

int_res <- rbindlist(lapply(names(outcomes), function(o) {
  r <- pool_contrasts(fits[[o]]$interaction, L_int)
  r[, outcome := o]
}))
int_res[, c("contrast", "level") := tstrsplit(term, "|", fixed = TRUE)]
int_res[contrast == "joint_m" & level == "Male:Normal", `:=`(se = NA, lcl = 0, ucl = 0, p = NA, df = NA, fmi = NA)]

# Global interaction tests and multiple-testing adjustment.
int_terms <- paste0("bmi_cat", bmi_nonref, ":sexFemale")
global <- rbindlist(lapply(names(outcomes), function(o) pool_wald(fits[[o]]$interaction, int_terms)[, outcome := o]))
global[, p_holm := p.adjust(p, "holm")]
int_res[contrast == "ratio", q_bh := p.adjust(p, "BH")]

# Flag categories resting on fewer than 30 unweighted cases in either sex.
int_res[, bmi := sub("^.*:", "", level)]
flag <- cases[, .(min_cases = min(n_cases)), by = .(outcome, bmi_cat)]
int_res[flag, on = .(outcome, bmi = bmi_cat), min_cases := i.min_cases]
int_res[, sparse := min_cases < 30]

fwrite(int_res[, .(outcome, contrast, level, log_or = est, se, or = round(exp(est), 3),
                   lcl = round(exp(lcl), 3), ucl = round(exp(ucl), 3), p, q_bh,
                   df = round(df, 1), fmi = round(fmi, 3), min_cases, sparse)],
       file.path(tab_dir, "models_interaction.csv"))
fwrite(global[, .(outcome, statistic = round(statistic, 3), df1, df2 = round(df2, 1), p, p_holm)],
       file.path(tab_dir, "interaction_tests.csv"))
fwrite(cases, file.path(tab_dir, "case_counts_sex_bmi.csv"))

# Table 3: adjusted ORs from the main-effects models -------------------------

term_label <- function(t) {
  lab <- c(bmi_cat = "BMI category", sex = "Sex", age_adj = "Age group", income = "Household income",
           educ = "Household education", smoking = "Smoking", alcohol = "Alcohol",
           phys_act = "Physical activity, MET-min/week")
  v <- names(lab)[vapply(names(lab), function(p) startsWith(t, p), logical(1))][1]
  list(group = lab[[v]], level = substring(t, nchar(v) + 1))
}
main_res[, c("group", "level") := rbindlist(lapply(term, term_label))]
main_res[, cell := paste0(fmt_or(est, lcl, ucl), fifelse(level == "Underweight", " &dagger;", ""))]
refs <- data.table(
  group = c("BMI category", "Sex", "Age group", "Household income", "Household education",
            "Smoking", "Alcohol", "Physical activity, MET-min/week"),
  level = c("Normal", "Male", levels(co$age_adj)[1], levels(co$income)[1], levels(co$educ)[1],
            levels(co$smoking)[1], levels(co$alcohol)[1], levels(co$phys_act)[1])
)
t3 <- dcast(main_res, group + level ~ outcome, value.var = "cell")
t3 <- rbind(t3, refs[, c(.SD, setNames(as.list(rep("1.00 (reference)", 4)), names(outcomes)))], fill = TRUE)
lev_order <- c("Normal", "Underweight", "Overweight", "Obese", "Male", "Female",
               levels(co$age_adj), levels(co$income), levels(co$educ), levels(co$smoking),
               levels(co$alcohol), levels(co$phys_act))
t3 <- t3[order(match(group, refs$group), match(level, lev_order))]
setcolorder(t3, c("group", "level", names(outcomes)))

t3_gt <- gt(t3, groupname_col = "group", rowname_col = "level") |>
  tab_header(title = "Table 3. Adjusted odds ratios from main-effects models",
             subtitle = sprintf("CCHS 2017-2018 PUMF, adults 18+ (n = %s). Odds ratio (95%% CI)", format(fits$cmm$main[[1]]$n, big.mark = ","))) |>
  cols_label(.list = as.list(outcomes)) |>
  fmt_markdown(columns = everything()) |>
  tab_source_note(paste(
    "Survey-weighted logistic regression, each outcome modelled separately with all variables shown.",
    sprintf("Bootstrap variance with 1,000 replicate weights; %d imputations pooled with Rubin's rules.", m),
    "BMI adjusted for self-report bias.",
    "&dagger; Fewer than 30 unweighted cases in at least one sex; estimate shown for completeness, not interpreted."
  ) |> html())
gtsave(t3_gt, file.path(tab_dir, "table3.html"))

# Table 4: sex x BMI interaction ----------------------------------------------

show <- int_res[contrast %in% c("men", "women", "ratio")]
show[, cell := paste0(fmt_or(est, lcl, ucl), fifelse(sparse, " &dagger;", ""))]
show[contrast == "ratio", cell := paste0(cell, "<br><small>p = ", fmt_p(p), "; q = ", fmt_p(q_bh), "</small>")]
t4 <- dcast(show, outcome + level ~ contrast, value.var = "cell")
t4 <- t4[order(match(outcome, names(outcomes)), match(level, bmi_nonref))]
t4[global, on = "outcome", global_p := paste0("Global interaction test: p = ", fmt_p(i.p), "; Holm-adjusted p = ", fmt_p(i.p_holm))]
t4[, outcome_label := paste0(outcomes[outcome], ": ", global_p)]
t4 <- t4[, .(outcome_label, level, men, women, ratio)]

t4_gt <- gt(t4, groupname_col = "outcome_label", rowname_col = "level") |>
  tab_header(title = "Table 4. Association of BMI category with each outcome, by sex",
             subtitle = "Odds ratio (95% CI) relative to normal weight in the same sex") |>
  cols_label(men = "Men", women = "Women", ratio = "Ratio of ORs (women / men)") |>
  fmt_markdown(columns = everything()) |>
  tab_source_note(paste(
    "From models with a sex x BMI category interaction, adjusted for age group, household income,",
    "household education, smoking, alcohol and physical activity. Ratio of ORs above 1: stronger association in women",
    "on the multiplicative scale. Global test: pooled Wald test (D1) of the three interaction terms; Holm adjustment",
    "across the four outcomes. q: Benjamini-Hochberg across the twelve ratio-of-OR tests.",
    "&dagger; Fewer than 30 unweighted cases in at least one sex; not interpreted.",
    "Joint odds ratios with a single reference group are in models_interaction.csv."
  ) |> html())
gtsave(t4_gt, file.path(tab_dir, "table4.html"))

print(global)
print(int_res[outcome == "cmm" & contrast %in% c("men", "women", "ratio"),
              .(contrast, level, or = fmt_or(est, lcl, ucl), p = fmt_p(p), q = fmt_p(q_bh), sparse)])

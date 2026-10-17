# Weighted descriptive tables on the observed (non-imputed) eligible cohort.
# Outputs: output/tables/table1.{csv,html}, output/tables/prevalence.{csv,html}

source("R/00_setup.R")
source("R/functions.R")
suppressPackageStartupMessages(library(gt))

co  <- readRDS(file.path(derived_dir, "cohort.rds"))
bsw <- readRDS(file.path(derived_dir, "bsw.rds"))
rw  <- rep_weights(co, bsw)
rm(bsw); invisible(gc())

labels <- c(
  bmi_cat = "BMI category (adjusted)", age_adj = "Age group", income = "Household income",
  educ = "Household education", smoking = "Smoking", alcohol = "Alcohol, past 12 months",
  phys_act = "Physical activity, MET-min/week", htn = "Hypertension", dm = "Diabetes",
  hd = "Heart disease", cmm = "Cardiometabolic multimorbidity"
)

# Conditions as factors so they appear in Table 1 with a missing row.
for (v in c("htn", "dm", "hd", "cmm")) co[, (paste0(v, "_f")) := factor(get(v), 0:1, c("No", "Yes"))]

# Table 1 ----------------------------------------------------------------------

t1_vars <- c("age_adj", "income", "educ", "smoking", "alcohol", "phys_act",
             "htn_f", "dm_f", "hd_f", "cmm_f")

co[, bmi_col := factor(bmi_cat)]
t1 <- rbind(
  rbindlist(lapply(c("bmi_cat", t1_vars), function(v) distribution(co, rw, v, "sex")))[, bmi_col := "All"],
  rbindlist(lapply(t1_vars, function(v) distribution(co, rw, v, c("sex", "bmi_col"))))
)
t1[, variable := sub("_f$", "", variable)]
t1[, cell := fifelse(flag %in% c("S", "F"), flag,
                     paste0(format(n_num, big.mark = ","), " (", fmt_pct(p), ")",
                            fifelse(flag == "E", " E", "")))]
# Drop rows for a "Missing" level that never occurs.
t1 <- t1[!(level == "Missing" & ave(n_num, variable, level, FUN = sum) == 0)]

# Column headers: unweighted n and weighted population (rounded to thousands).
hdr <- rbind(
  co[, .(n = .N, pop = sum(WTS_M)), by = sex][, bmi_col := "All"],
  co[!is.na(bmi_col), .(n = .N, pop = sum(WTS_M)), by = .(sex, bmi_col)]
)
hdr[, label := sprintf("n = %s; weighted N = %s", format(n, big.mark = ","),
                       format(round(pop, -3), big.mark = ",", scientific = FALSE))]

col_order <- c("All", levels(co$bmi_cat))
t1[, column := paste(sex, bmi_col, sep = "_")]
t1_wide <- dcast(t1, variable + level ~ column, value.var = "cell", fill = "")
cols <- c(paste0("Male_", col_order), paste0("Female_", col_order))
var_order <- c("bmi_cat", sub("_f$", "", t1_vars))
lev_order <- unique(t1[, .(variable, level)])
t1_wide <- t1_wide[lev_order, on = .(variable, level)]
t1_wide <- t1_wide[order(match(variable, var_order))]
setcolorder(t1_wide, c("variable", "level", cols))
t1_wide[, variable := labels[variable]]

fwrite(t1[, .(variable, level, sex, bmi_col, n = n_num, n_domain = n_den,
              weighted_pct = round(100 * p, 1), se_pct = round(100 * se, 2),
              cv = round(cv, 3), flag)],
       file.path(tab_dir, "table1.csv"))

spanner_label <- function(s) {
  h <- hdr[sex == s]
  setNames(h$label, paste(h$sex, h$bmi_col, sep = "_"))
}
col_labels <- c(spanner_label("Male"), spanner_label("Female"))
col_labels <- setNames(
  sprintf("%s<br><small>%s</small>", sub("^.*_", "", names(col_labels)), col_labels),
  names(col_labels)
)
col_labels <- sub("^All<br>", "Total<br>", col_labels)

t1_gt <- gt(t1_wide, groupname_col = "variable", rowname_col = "level") |>
  tab_header(title = "Table 1. Characteristics of adults by sex and BMI category",
             subtitle = "CCHS 2017-2018 PUMF. Unweighted n (weighted %)") |>
  tab_spanner("Men", columns = starts_with("Male_")) |>
  tab_spanner("Women", columns = starts_with("Female_")) |>
  cols_label(.list = lapply(col_labels, html)) |>
  tab_source_note(paste(
    "Percentages are weighted with WTS_M; columns sum to 100% within each variable.",
    "BMI adjusted for self-report bias. Respondents with missing BMI appear in the Total columns only.",
    "E: coefficient of variation 15.1-35.0%, interpret with caution.",
    "S: fewer than 30 respondents, not released. F: coefficient of variation above 35.0%, not released."
  ))
gtsave(t1_gt, file.path(tab_dir, "table1.html"))

# Weighted prevalence ----------------------------------------------------------

conds <- c("htn", "dm", "hd", "cmm")
domains <- list(
  overall       = NULL,
  sex           = "sex",
  sex_bmi       = c("sex", "bmi_cat"),
  sex_age       = c("sex", "age_grp"),
  sex_age_bmi   = c("sex", "age_grp", "bmi_cat")
)

prev <- rbindlist(lapply(names(domains), function(nm) {
  rbindlist(lapply(conds, function(v) prevalence(co, rw, v, domains[[nm]])))[, domain := nm]
}), fill = TRUE)
for (v in c("sex", "bmi_cat", "age_grp")) prev[is.na(get(v)), (v) := "All"]
prev[, `:=`(sex = factor(sex, c("All", "Male", "Female")),
            age_grp = factor(age_grp, c("All", levels(co$age_grp))),
            bmi_cat = factor(bmi_cat, c("All", levels(co$bmi_cat))))]
setorder(prev, domain, sex, age_grp, bmi_cat)

fwrite(prev[, .(domain, sex, age_grp, bmi_cat, condition = variable,
                n_cases = n_num, n = n_den,
                prevalence_pct = round(100 * p, 1), lcl_pct = round(100 * lcl, 1),
                ucl_pct = round(100 * ucl, 1), se_pct = round(100 * se, 2),
                cv = round(cv, 3), flag)],
       file.path(tab_dir, "prevalence.csv"))

prev[, cell := fmt_prev(p, lcl, ucl, flag)]
prev_wide <- dcast(prev[domain != "sex_age_bmi"], domain + sex + age_grp + bmi_cat ~ variable, value.var = "cell")
prev_wide[, domain := factor(domain, names(domains))]
setorder(prev_wide, domain, sex, age_grp, bmi_cat)
prev_wide[, group := fcase(
  domain == "overall", "All adults",
  domain == "sex",     "By sex",
  domain == "sex_bmi", paste0(sex, ", by BMI category"),
  domain == "sex_age", paste0(sex, ", by age group")
)]
prev_wide[, row := fcase(
  domain == "overall", "All",
  domain == "sex",     as.character(sex),
  domain == "sex_bmi", as.character(bmi_cat),
  domain == "sex_age", as.character(age_grp)
)]

prev_gt <- gt(prev_wide[, .(group, row, htn, dm, hd, cmm)], groupname_col = "group", rowname_col = "row") |>
  tab_header(title = "Table 2. Weighted prevalence of cardiometabolic conditions",
             subtitle = "CCHS 2017-2018 PUMF, adults 18+. Percent (95% CI)") |>
  cols_label(htn = labels[["htn"]], dm = labels[["dm"]], hd = labels[["hd"]], cmm = labels[["cmm"]]) |>
  tab_source_note(paste(
    "Bootstrap variance with 1,000 replicate weights; 95% CI on the logit scale.",
    "Multimorbidity: two or more of hypertension, diabetes and heart disease.",
    "Diabetes excludes diabetes diagnosed only during pregnancy. BMI adjusted for self-report bias.",
    "E: coefficient of variation 15.1-35.0%. S: fewer than 30 respondents, not released.",
    "F: coefficient of variation above 35.0%, not released.",
    "Prevalence by sex, age group and BMI category together is in prevalence.csv."
  ))
gtsave(prev_gt, file.path(tab_dir, "prevalence.html"))

print(prev[domain %in% c("overall", "sex", "sex_bmi") & variable == "cmm",
           .(sex, bmi_cat, n_cases = n_num, n = n_den, cell)])
print(prev[, .N, by = flag])
print(t1[, .N, by = flag])

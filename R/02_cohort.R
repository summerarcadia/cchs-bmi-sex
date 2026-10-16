# Define the study cohort, recode analysis variables and derive the outcome.
# Codes follow the 2017-2018 PUMF Data Dictionary and Derived Variables documents.
# Outputs: data/derived/cohort.rds (not tracked), output/tables/flow_counts.csv,
#          output/tables/missing_by_variable.csv

source("R/00_setup.R")

d <- readRDS(file.path(derived_dir, "pumf_selected.rds"))

yes_no <- function(x) fifelse(x == 1L, 1L, fifelse(x == 2L, 0L, NA_integer_))
na_codes <- function(x, codes) replace(x, x %in% codes, NA)

# Population flags -----------------------------------------------------------

d[, adult := DHHGAGE >= 3L]                       # codes 03-16 = 18 and older
# HWTDGBCC = 6 (valid skip) for an adult means MAC_025 = 1, i.e. currently
# pregnant. MAC_025 itself is not on the PUMF.
d[, pregnant := adult & HWTDGBCC == 6L]

# Demographics ---------------------------------------------------------------

d[, sex := factor(DHH_SEX, 1:2, c("Male", "Female"))]

d[, age_grp := fcase(
  DHHGAGE %between% c(3L, 8L),   "18-44",
  DHHGAGE %between% c(9L, 12L),  "45-64",
  DHHGAGE %between% c(13L, 16L), "65+"
)]
d[, age_grp := factor(age_grp, c("18-44", "45-64", "65+"))]

# Finer age for adjustment, built from the 5-year PUMF groups.
d[, age_adj := fcase(
  DHHGAGE %between% c(3L, 6L),   "18-34",
  DHHGAGE %between% c(7L, 8L),   "35-44",
  DHHGAGE %between% c(9L, 10L),  "45-54",
  DHHGAGE %between% c(11L, 12L), "55-64",
  DHHGAGE %between% c(13L, 14L), "65-74",
  DHHGAGE %between% c(15L, 16L), "75+"
)]
d[, age_adj := factor(age_adj, c("18-34", "35-44", "45-54", "55-64", "65-74", "75+"))]

d[, province := factor(GEO_PRV)]
d[, proxy := ADM_PRX == 1L]

# BMI ------------------------------------------------------------------------
# Primary exposure: BMI adjusted for self-report bias (HWTDGBCC, from HWTDGCOR,
# sex-specific correction equations of Connor Gorber et al. 2008).
# Self-reported BMI (HWTDGISW / HWTDGBMI) is kept for sensitivity analysis.

bmi_levels <- c("Underweight", "Normal", "Overweight", "Obese")
d[, bmi_cat    := factor(na_codes(HWTDGBCC, c(6L, 9L)), 1:4, bmi_levels)]
d[, bmi_cat_sr := factor(na_codes(HWTDGISW, c(6L, 9L)), 1:4, bmi_levels)]
d[, bmi_adj    := na_codes(HWTDGCOR, c(999.96, 999.99))]
d[, bmi_sr     := na_codes(HWTDGBMI, c(999.96, 999.99))]

# Conditions -----------------------------------------------------------------

d[, htn := yes_no(CCC_065)]
d[, hd  := yes_no(CCC_085)]

# Diabetes, excluding diabetes diagnosed only during pregnancy
# (CCC_105 = 1 pregnant at first diagnosis, CCC_110 = 2 never told otherwise).
d[, dm := yes_no(CCC_095)]
d[dm == 1L & CCC_105 == 1L & CCC_110 == 2L, dm := 0L]
d[dm == 1L & CCC_105 == 1L & CCC_110 %in% c(7L, 8L, 9L), dm := NA_integer_]

# Cardiometabolic multimorbidity: two or more of hypertension, diabetes,
# heart disease. Determined when the known conditions settle it.
d[, n_cond := rowSums(cbind(htn, dm, hd) == 1L, na.rm = TRUE)]
d[, n_miss := rowSums(is.na(cbind(htn, dm, hd)))]
d[, cmm := fcase(
  n_cond >= 2L,          1L,
  n_cond + n_miss < 2L,  0L,
  default = NA_integer_
)]

# Covariates -----------------------------------------------------------------

d[, educ := factor(na_codes(EHG2DVH3, 9L), 1:3,
                   c("Less than secondary", "Secondary", "Post-secondary"))]

d[, income := factor(na_codes(INCDGHH, 9L), 1:5,
                     c("<20k", "20-39k", "40-59k", "60-79k", "80k+"))]

d[, smoking := fcase(
  SMKDVSTY %in% c(1L, 2L), "Current",
  SMKDVSTY %in% c(3L, 4L), "Former",
  SMKDVSTY %in% c(5L, 6L), "Never"
)]
d[, smoking := factor(smoking, c("Never", "Former", "Current"))]

d[, alcohol := factor(na_codes(ALCDVTTM, 9L), c(3L, 2L, 1L),
                      c("None in past 12 mo", "Occasional", "Regular"))]

# PAADVWHO: weekly activity volume in MET-minutes (Derived Variables, p. 116).
d[, phys_act := factor(na_codes(PAADVWHO, c(6L, 9L)), 1:4,
                       c(">=900", "450-899", "1-449", "0"))]
d[, phys_act := factor(phys_act, c("0", "1-449", "450-899", ">=900"))]

# Flow -----------------------------------------------------------------------

n_all      <- nrow(d)
n_youth    <- d[adult == FALSE, .N]
n_preg     <- d[pregnant == TRUE, .N]
cohort     <- d[adult == TRUE & pregnant == FALSE]
n_cohort   <- nrow(cohort)

model_vars <- c("cmm", "bmi_cat", "sex", "age_adj", "income", "educ", "smoking", "alcohol", "phys_act")
cohort[, complete := complete.cases(.SD), .SDcols = model_vars]

flow <- data.table(
  step = c("All PUMF respondents",
           "Excluded: aged 12-17",
           "Excluded: currently pregnant",
           "Eligible cohort (adults 18+)",
           "  of which missing outcome",
           "  of which missing BMI category (outcome known)",
           "  of which missing covariate only",
           "Complete-case sample"),
  n = c(n_all, n_youth, n_preg, n_cohort,
        cohort[is.na(cmm), .N],
        cohort[!is.na(cmm) & is.na(bmi_cat), .N],
        cohort[!is.na(cmm) & !is.na(bmi_cat) & !complete, .N],
        cohort[complete == TRUE, .N])
)
stopifnot(n_all - n_youth - n_preg == n_cohort)
fwrite(flow, file.path(tab_dir, "flow_counts.csv"))

# Unweighted missingness in the eligible cohort, by sex.
miss_vars <- c("htn", "dm", "hd", "cmm", "bmi_cat", "bmi_cat_sr", "bmi_adj",
               "income", "educ", "smoking", "alcohol", "phys_act")
miss <- rbindlist(lapply(miss_vars, function(v) {
  cohort[, .(variable = v, n_missing = sum(is.na(get(v))), n = .N), by = sex]
}))
miss[, pct_missing := round(100 * n_missing / n, 2)]
miss <- dcast(miss, variable ~ sex, value.var = c("n_missing", "pct_missing"))
setcolorder(miss, c("variable", "n_missing_Male", "pct_missing_Male", "n_missing_Female", "pct_missing_Female"))
miss <- miss[match(miss_vars, variable)]
fwrite(miss, file.path(tab_dir, "missing_by_variable.csv"))

# Proxy interviews skip the height and weight module, so BMI is structurally
# missing for them.
proxy_bmi <- cohort[, .(n = .N, bmi_missing = sum(is.na(bmi_cat))), by = proxy]

keep_out <- c("ADM_RNO", "WTS_M", "sex", "age_grp", "age_adj", "province", "proxy",
              "bmi_cat", "bmi_cat_sr", "bmi_adj", "bmi_sr",
              "htn", "dm", "hd", "n_cond", "cmm",
              "educ", "income", "smoking", "alcohol", "phys_act", "complete")
saveRDS(cohort[, ..keep_out], file.path(derived_dir, "cohort.rds"))

print(flow)
print(miss)
print(proxy_bmi)
print(cohort[, .N, keyby = .(sex, bmi_cat)])
print(cohort[, .N, keyby = .(age_grp, n_cond)])

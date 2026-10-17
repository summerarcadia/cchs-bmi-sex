# E-values for unmeasured confounding (VanderWeele and Ding 2017) for the
# sex-specific odds ratios of overweight and obesity against normal weight
# (R/07_model_tables.R). An outcome is treated as rare when its weighted
# prevalence in the cohort is below 15% (Table 2); otherwise the odds ratio is
# converted to an approximate risk ratio by the square-root transformation.
#
# Output: output/tables/evalues.csv

source("R/00_setup.R")
suppressPackageStartupMessages(library(EValue))

ors  <- fread(file.path(tab_dir, "models_interaction.csv"))
prev <- fread(file.path(tab_dir, "prevalence.csv"))[domain == "overall"]
rare <- setNames(prev$prevalence_pct < 15, prev$condition)

ev <- ors[contrast %in% c("men", "women") & level %in% c("Overweight", "Obese")]
ev[, c("evalue", "evalue_ci") := {
  r <- t(mapply(function(o, l, u, out) {
    e <- evalues.OR(o, l, u, rare = rare[[out]])
    ci <- e["E-values", c("lower", "upper")]
    c(e["E-values", "point"], ci[!is.na(ci)][1])
  }, or, lcl, ucl, outcome))
  list(r[, 1], r[, 2])
}]
ev[, rare_outcome := rare[outcome]]

fwrite(ev[, .(outcome, sex = contrast, bmi = level, or, lcl, ucl, rare_outcome,
              evalue = round(evalue, 2), evalue_ci = round(evalue_ci, 2))],
       file.path(tab_dir, "evalues.csv"))
print(ev[, .(outcome, sex = contrast, bmi = level, or, lcl, ucl, rare_outcome, evalue = round(evalue, 2), evalue_ci = round(evalue_ci, 2))])

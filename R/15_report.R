# Build the manuscript and one-page summary as Word documents.
#
# The text lives in manuscript/manuscript.md and manuscript/summary.md. Each
# {{key}} is replaced with a value computed here from the files in
# output/tables, so every number in the documents comes from the analysis
# outputs. {{table:name}} and {{figure:file}} lines insert a table or figure.
# An unknown key stops the build.
#
# Outputs: output/report/manuscript.docx, output/report/summary.docx

source("R/00_setup.R")
source("R/functions.R")
suppressPackageStartupMessages({
  library(officer)
  library(flextable)
})

report_dir <- here("output", "report")
dir.create(report_dir, showWarnings = FALSE)
rd <- function(f) fread(file.path(tab_dir, f))

# Formatting --------------------------------------------------------------------

n_fmt <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
d1 <- function(x) formatC(x, format = "f", digits = 1)
d2 <- function(x) formatC(x, format = "f", digits = 2)
ci2 <- function(e, l, u) sprintf("%s [95%% CI %s to %s]", d2(e), d2(l), d2(u))
ci1 <- function(e, l, u) sprintf("%s [95%% CI %s to %s]", d1(e), d1(l), d1(u))
p_txt <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %s", formatC(p, format = "f", digits = 3))

# Inputs ------------------------------------------------------------------------

flow  <- rd("flow_counts.csv"); flow[, step := trimws(step)]
fl    <- setNames(flow$n, flow$step)
miss  <- rd("missing_by_variable.csv")
prev  <- rd("prevalence.csv")
ors   <- rd("models_interaction.csv")
tests <- rd("interaction_tests.csv")
main  <- rd("models_main.csv")
add   <- rd("additive_interaction.csv")
std   <- rd("standardized_prevalence.csv")
ast   <- rd("age_strata_additive.csv")
atest <- rd("age_strata_tests.csv")
sens  <- rd("sensitivity.csv")
ps    <- rd("propensity.csv")
bal   <- rd("propensity_balance.csv")
ev    <- rd("evalues.csv")
co    <- readRDS(file.path(derived_dir, "cohort.rds"))
fits_m <- length(readRDS(file.path(derived_dir, "imputed.rds")))

pv <- function(dom, cond, s = "All", b = "All", a = "All") {
  prev[domain == dom & condition == cond & sex == s & bmi_cat == b & age_grp == a]
}
# Argument names are prefixed so they cannot collide with column names.
orv <- function(.o, .c, .l) ors[outcome == .o & contrast == .c & level == .l]
addv <- function(.o, .b, .m) add[outcome == .o & bmi == .b & measure == .m]
stdv <- function(.o, .b, .m) std[outcome == .o & bmi == .b & measure == .m]
astv <- function(.o, .s, .b, .m) ast[outcome == .o & stratum == .s & bmi == .b & measure == .m]
sensv <- function(.v, .b, .m) sens[variant == .v & bmi == .b & measure == .m]
psv <- function(.a, .o, .b, .m) ps[approach == .a & outcome == .o & bmi == .b & measure == .m]

v <- list()
v$n_all <- n_fmt(fl[["All PUMF respondents"]])
v$n_youth <- n_fmt(fl[["Excluded: aged 12-17"]])
v$n_preg <- n_fmt(fl[["Excluded: currently pregnant"]])
v$n_cohort <- n_fmt(fl[["Eligible cohort (adults 18+)"]])
v$n_complete <- n_fmt(fl[["Complete-case sample"]])
v$pct_complete <- d1(100 * fl[["Complete-case sample"]] / fl[["Eligible cohort (adults 18+)"]])
v$m <- fits_m
v$miss_bmi_men <- d1(miss[variable == "bmi_cat", pct_missing_Male])
v$miss_bmi_women <- d1(miss[variable == "bmi_cat", pct_missing_Female])
v$n_proxy <- n_fmt(co[, sum(proxy)])

# Prevalence.
x <- pv("overall", "cmm"); v$prev_cmm <- sprintf("%s%% [95%% CI %s%% to %s%%]", d1(x$prevalence_pct), d1(x$lcl_pct), d1(x$ucl_pct))
for (s in c("Male", "Female")) {
  key <- tolower(substr(s, 1, 1))
  x <- pv("sex", "cmm", s); v[[paste0("prev_cmm_", key)]] <- d1(x$prevalence_pct)
  for (b in c("Normal", "Obese")) {
    x <- pv("sex_bmi", "cmm", s, b); v[[paste0("prev_cmm_", key, "_", tolower(b))]] <- d1(x$prevalence_pct)
  }
}
x <- pv("sex_age", "cmm", "Male", "All", "65+"); v$prev_cmm_m_65 <- d1(x$prevalence_pct)
x <- pv("sex_age", "cmm", "Female", "All", "65+"); v$prev_cmm_f_65 <- d1(x$prevalence_pct)
v$n_uw_cases_m <- ors[outcome == "cmm" & contrast == "men" & level == "Underweight", min_cases]

# Odds ratios.
for (o in c("cmm", "dm", "htn", "hd")) for (b in c("Overweight", "Obese")) {
  for (cc in c("men", "women", "ratio")) {
    x <- orv(o, cc, b); v[[paste("or", o, cc, tolower(b), sep = "_")]] <- ci2(x$or, x$lcl, x$ucl)
  }
}
for (o in c("cmm", "dm", "htn", "hd")) {
  x <- tests[outcome == o]
  v[[paste0("pint_", o)]] <- p_txt(x$p); v[[paste0("pholm_", o)]] <- p_txt(x$p_holm)
}
x <- main[outcome == "cmm" & term == "sexFemale"]; v$or_female_main <- ci2(x$or, x$lcl, x$ucl)
x <- orv("cmm", "sex", "Normal"); v$or_female_normal <- ci2(x$or, x$lcl, x$ucl)

# Additive interaction, whole cohort.
for (o in c("cmm", "dm", "htn", "hd")) for (b in c("Overweight", "Obese")) {
  bb <- tolower(b)
  for (mm in c("reri", "ap", "si")) {
    x <- addv(o, b, mm)
    v[[paste(mm, o, bb, sep = "_")]] <- if (is.na(x$lcl)) paste0(d2(x$value), " [CI not estimable]") else ci2(x$value, x$lcl, x$ucl)
  }
  for (mm in c("rd_female", "rd_male", "rd_diff")) {
    x <- stdv(o, b, mm); v[[paste(mm, o, bb, sep = "_")]] <- ci1(x$pct, x$lcl_pct, x$ucl_pct)
    v[[paste(mm, o, bb, "pt", sep = "_")]] <- d1(x$pct)
  }
}

# Age strata (multimorbidity, obesity).
for (s in c("18-44", "45-64", "65+")) {
  key <- gsub("[^0-9]", "", s)
  for (mm in c("rd_female", "rd_male")) {
    x <- astv("cmm", s, "Obese", mm); v[[paste(mm, "obese", key, sep = "_")]] <- d1(100 * x$value)
  }
  x <- astv("cmm", s, "Obese", "rd_diff"); v[[paste("rd_diff_obese", key, sep = "_")]] <- ci1(100 * x$value, 100 * x$lcl, 100 * x$ucl)
  x <- astv("cmm", s, "Obese", "reri"); v[[paste("reri_obese", key, sep = "_")]] <- ci2(x$value, x$lcl, x$ucl)
}
ref_young <- co[age_grp == "18-44" & bmi_cat == "Normal" & !is.na(cmm), .(cases = sum(cmm)), by = sex]
v$n_ref_cases_young_m <- ref_young[sex == "Male", cases]; v$n_ref_cases_young_f <- ref_young[sex == "Female", cases]
x <- atest[outcome == "dm"]; v$pint_dm_45 <- p_txt(x[stratum == "45-64", p]); v$pint_dm_65 <- p_txt(x[stratum == "65+", p])

# Sensitivity.
for (vv in c("complete_case", "excl_underweight", "self_report")) {
  x <- sensv(vv, "Obese", "reri"); v[[paste0("reri_", vv)]] <- ci2(x$value, x$lcl, x$ucl)
  x <- sensv(vv, "Obese", "rd_diff"); v[[paste0("rd_diff_", vv)]] <- ci1(100 * x$value, 100 * x$lcl, 100 * x$ucl)
}
x <- sensv("obesity_classes", "ObeseI", "reri"); v$reri_class1 <- ci2(x$value, x$lcl, x$ucl)
x <- sensv("obesity_classes", "ObeseII", "reri"); v$reri_class2 <- ci2(x$value, x$lcl, x$ucl)
x <- sensv("obesity_classes", "ObeseII", "rd_diff"); v$rd_diff_class2 <- ci1(100 * x$value, 100 * x$lcl, 100 * x$ucl)
x <- psv("doubly_robust", "cmm", "Obese", "reri"); v$reri_ps_dr <- ci2(x$value, x$lcl, x$ucl)
x <- psv("doubly_robust", "cmm", "Obese", "rd_diff"); v$rd_diff_ps_dr <- ci1(100 * x$value, 100 * x$lcl, 100 * x$ucl)
x <- psv("marginal", "cmm", "Obese", "reri"); v$reri_ps_marg <- ci2(x$value, x$lcl, x$ucl)
x <- psv("marginal", "cmm", "Obese", "rd_diff"); v$rd_diff_ps_marg <- ci1(100 * x$value, 100 * x$lcl, 100 * x$ucl)
v$smd_unadj <- d2(max(bal$max_smd_unadjusted)); v$smd_wtd <- d2(max(bal$max_smd_weighted))

# E-values.
x <- ev[outcome == "cmm" & bmi == "Obese"]
v$ev_cmm_obese_m <- d2(x[sex == "men", evalue]); v$ev_cmm_obese_m_ci <- d2(x[sex == "men", evalue_ci])
v$ev_cmm_obese_f <- d2(x[sex == "women", evalue]); v$ev_cmm_obese_f_ci <- d2(x[sex == "women", evalue_ci])
v$ev_min <- d2(min(ev$evalue_ci))
v$ev_min_row <- with(ev[which.min(evalue_ci)], sprintf("%s, %s, %s", c(cmm = "multimorbidity", htn = "hypertension", dm = "diabetes", hd = "heart disease")[outcome], sex, tolower(bmi)))

bad <- names(v)[vapply(v, function(x) length(x) != 1L || is.na(x), logical(1))]
if (length(bad)) stop("Values missing or not scalar: ", paste(bad, collapse = ", "))

# Tables -------------------------------------------------------------------------

std_ft <- function(dt, header_map, widths = NULL) {
  ft <- flextable(dt) |> set_header_labels(values = header_map) |>
    fontsize(size = 8, part = "all") |> font(fontname = "Arial", part = "all") |>
    padding(padding = 1.5, part = "all") |> bold(part = "header") |>
    border_remove() |> hline_top(border = fp_border(width = 1), part = "header") |>
    hline_bottom(border = fp_border(width = 0.75), part = "header") |>
    hline_bottom(border = fp_border(width = 1), part = "body") |>
    valign(valign = "top", part = "all")
  if (!is.null(widths)) ft <- width(ft, width = widths)
  ft
}

# Table 1: characteristics by sex (Total columns of output/tables/table1.csv).
t1 <- rd("table1.csv")[bmi_col == "All"]
t1_vars <- c(bmi_cat = "BMI category (adjusted)", age_adj = "Age group", income = "Household income",
             educ = "Household education", smoking = "Smoking", alcohol = "Alcohol, past 12 months",
             phys_act = "Physical activity, MET-min/week", htn = "Hypertension", dm = "Diabetes",
             hd = "Heart disease", cmm = "Multimorbidity")
t1 <- t1[variable %in% names(t1_vars) & !(variable %in% c("htn", "dm", "hd", "cmm") & level %in% c("No", "Missing"))]
t1[, cell := fifelse(flag %in% c("S", "F"), flag, paste0(n_fmt(n), " (", d1(weighted_pct), ")", fifelse(flag == "E", " E", "")))]
t1w <- dcast(t1, variable + level ~ sex, value.var = "cell")
lev_order <- unique(rd("table1.csv")[, .(variable, level)])[, ord := .I]
t1w <- merge(t1w, lev_order, by = c("variable", "level"))[order(match(variable, names(t1_vars)), ord)]
t1w[, variable := fifelse(duplicated(variable), "", t1_vars[variable])]
t1w <- t1w[, .(variable, level, Female, Male)]
hdr_n <- rd("table1.csv")[bmi_col == "All" & variable == "age_adj", .(n = sum(n)), by = sex]
ft1 <- std_ft(t1w, c(variable = "Characteristic", level = "",
                     Female = sprintf("Women (n = %s)", n_fmt(hdr_n[sex == "Female", n])),
                     Male = sprintf("Men (n = %s)", n_fmt(hdr_n[sex == "Male", n]))),
              widths = c(1.9, 1.6, 1.4, 1.4))

# Table 2: sex-specific odds ratios and ratio of odds ratios.
oc <- c(cmm = "Multimorbidity", htn = "Hypertension", dm = "Diabetes", hd = "Heart disease")
t2 <- ors[contrast %in% c("men", "women", "ratio") & level %in% c("Overweight", "Obese")]
t2[, cell := sprintf("%s (%s-%s)", d2(or), d2(lcl), d2(ucl))]
t2w <- dcast(t2, outcome + level ~ contrast, value.var = "cell")
t2w[tests, on = "outcome", p_int := fifelse(i.p < 0.001, "<0.001", formatC(i.p, format = "f", digits = 3))]
t2w[tests, on = "outcome", p_holm := fifelse(i.p_holm < 0.001, "<0.001", formatC(i.p_holm, format = "f", digits = 3))]
t2w <- t2w[order(match(outcome, names(oc)), match(level, c("Overweight", "Obese")))]
t2w[level == "Obese", `:=`(p_int = "", p_holm = "")]
t2w[, outcome := fifelse(duplicated(outcome), "", oc[outcome])]
t2w <- t2w[, .(outcome, level, women, men, ratio, p_int, p_holm)]
ft2 <- std_ft(t2w, c(outcome = "Outcome", level = "BMI category", women = "OR, women", men = "OR, men",
                     ratio = "Ratio of ORs (women/men)", p_int = "Interaction p", p_holm = "Holm p"),
              widths = c(1.1, 0.85, 1.05, 1.05, 1.15, 0.7, 0.6))

# Table 3: additive interaction.
t3 <- merge(
  dcast(std[measure %in% c("rd_female", "rd_male", "rd_diff")][, cell := sprintf("%s (%s to %s)", d1(pct), d1(lcl_pct), d1(ucl_pct))],
        outcome + bmi ~ measure, value.var = "cell"),
  dcast(add[measure %in% c("reri", "ap")][, cell := sprintf("%s (%s to %s)", d2(value), d2(lcl), d2(ucl))],
        outcome + bmi ~ measure, value.var = "cell"),
  by = c("outcome", "bmi"))
t3 <- t3[order(match(outcome, names(oc)), match(bmi, c("Overweight", "Obese")))]
t3[, outcome := fifelse(duplicated(outcome), "", oc[outcome])]
t3 <- t3[, .(outcome, bmi, rd_female, rd_male, rd_diff, reri, ap)]
ft3 <- std_ft(t3, c(outcome = "Outcome", bmi = "BMI category", rd_female = "Women", rd_male = "Men",
                    rd_diff = "Women minus men", reri = "RERI", ap = "AP"),
              widths = c(1.05, 0.85, 1.0, 1.0, 1.1, 1.05, 1.0)) |>
  add_header_row(values = c("", "", "Standardized prevalence difference (percentage points)", "Odds ratio scale"),
                 colwidths = c(1, 1, 3, 2)) |>
  fontsize(size = 8, part = "header") |> bold(part = "header") |> align(i = 1, align = "center", part = "header")

tables <- list(table1 = ft1, table2 = ft2, table3 = ft3)

# Renderer -----------------------------------------------------------------------

fill <- function(line) {
  keys <- regmatches(line, gregexpr("\\{\\{[a-z0-9_]+\\}\\}", line))[[1]]
  for (k in keys) {
    nm <- gsub("[{}]", "", k)
    if (is.null(v[[nm]])) stop("Unknown key in template: ", nm)
    line <- sub(k, as.character(v[[nm]]), line, fixed = TRUE)
  }
  line
}

render <- function(template, out, title) {
  lines <- readLines(template, warn = FALSE)
  stopifnot(!any(grepl("[^ -~]", lines)))
  doc <- read_docx()
  doc <- set_doc_properties(doc, title = title, creator = "Hayley So", description = "",
                            subject = "", keywords = "", created = Sys.time())
  body_style <- fp_text(font.size = 10.5, font.family = "Arial")
  head_style <- list(
    h1 = fp_text(font.size = 15, font.family = "Arial", bold = TRUE),
    h2 = fp_text(font.size = 12, font.family = "Arial", bold = TRUE),
    h3 = fp_text(font.size = 10.5, font.family = "Arial", bold = TRUE, italic = TRUE)
  )
  add_par <- function(doc, txt, style = "Normal") {
    body_add_fpar(doc, fpar(ftext(txt, body_style), fp_p = fp_par(line_spacing = 1.15, padding.bottom = 4)), style = style)
  }
  add_head <- function(doc, txt, level) {
    body_add_fpar(doc, fpar(ftext(txt, head_style[[level]]),
                            fp_p = fp_par(padding.top = if (level == "h1") 0 else 8, padding.bottom = 3)))
  }
  bullet <- intToUtf8(8226)
  add_bullet <- function(doc, txt) {
    body_add_fpar(doc, fpar(ftext(paste0(bullet, "  ", txt), body_style),
                            fp_p = fp_par(line_spacing = 1.15, padding.bottom = 2, padding.left = 14)))
  }
  buf <- character()
  flush <- function(doc) {
    if (length(buf)) { doc <- add_par(doc, paste(buf, collapse = " ")); buf <<- character() }
    doc
  }
  for (ln in lines) {
    ln <- fill(ln)
    if (grepl("^\\s*$", ln)) { doc <- flush(doc); next }
    if (startsWith(ln, "# ")) { doc <- flush(doc); doc <- add_head(doc, substring(ln, 3), "h1"); next }
    if (startsWith(ln, "## ")) { doc <- flush(doc); doc <- add_head(doc, substring(ln, 4), "h2"); next }
    if (startsWith(ln, "### ")) { doc <- flush(doc); doc <- add_head(doc, substring(ln, 5), "h3"); next }
    if (grepl("^\\[\\[table:", ln)) {
      doc <- flush(doc); nm <- sub("^\\[\\[table:([a-z0-9_]+)\\]\\]$", "\\1", ln)
      doc <- body_add_flextable(doc, tables[[nm]], align = "left"); next
    }
    if (grepl("^\\[\\[figure:", ln)) {
      doc <- flush(doc); spec <- sub("^\\[\\[figure:(.+)\\]\\]$", "\\1", ln)
      parts <- strsplit(spec, "|", fixed = TRUE)[[1]]
      doc <- body_add_img(doc, src = file.path(fig_dir, parts[1]), width = as.numeric(parts[2]), height = as.numeric(parts[3])); next
    }
    if (ln == "[[pagebreak]]") { doc <- flush(doc); doc <- body_add_break(doc); next }
    if (startsWith(ln, "- ")) { doc <- flush(doc); doc <- add_bullet(doc, substring(ln, 3)); next }
    if (grepl("^[0-9]+\\. ", ln)) { doc <- flush(doc); doc <- add_par(doc, ln); next }
    buf <- c(buf, ln)
  }
  doc <- flush(doc)
  print(doc, target = out)
}

render(here("manuscript", "manuscript.md"), file.path(report_dir, "manuscript.docx"),
       "Sex differences in the association between BMI and cardiometabolic multimorbidity in Canadian adults")
render(here("manuscript", "summary.md"), file.path(report_dir, "summary.docx"),
       "Summary: BMI, sex and cardiometabolic multimorbidity in Canadian adults")
cat("Documents written to", report_dir, "\n")

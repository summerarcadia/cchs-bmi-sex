# Multiple imputation of missing exposure, outcome components and covariates.
#
# Imputation is run separately within each sex x age group stratum so that
# associations that differ by sex and age (including sex x BMI) are carried
# into the imputed values. Continuous adjusted BMI is imputed by predictive
# mean matching; BMI categories, self-reported BMI and multimorbidity are
# derived from the completed data.
#
# Outputs (not tracked): data/derived/imputed.rds (list of m data.tables),
#                        data/derived/mids_by_stratum.rds
# Tracked: output/figures/mi_trace.png, output/tables/mi_observed_vs_imputed.csv,
#          output/tables/prevalence_cmm_mi.csv

source("R/00_setup.R")
source("R/functions.R")
suppressPackageStartupMessages({
  library(mice)
  library(parallel)
  library(ggplot2)
})

m      <- as.integer(Sys.getenv("MI_M", "20"))
maxit  <- as.integer(Sys.getenv("MI_MAXIT", "10"))
n_core <- max(1L, detectCores() - 2L)
seed   <- 20261016L

co <- readRDS(file.path(derived_dir, "cohort.rds"))

# Self-reported BMI categories must be recoverable from adjusted BMI by
# inverting the sex-specific correction (Derived Variables, HWTDGCOR).
corr <- data.table(sex = factor(c("Male", "Female"), levels(co$sex)),
                   a = c(-1.07575, -0.12374), b = c(1.07592, 1.05129))
bmi_breaks <- c(-Inf, 18.5, 25, 30, Inf)
bmi_levels <- levels(co$bmi_cat)
derive_bmi <- function(dt) {
  dt[corr, on = "sex", bmi_sr_inv := (bmi_adj - i.a) / i.b]
  dt[, bmi_cat_d    := cut(bmi_adj,    bmi_breaks, bmi_levels, right = FALSE)]
  dt[, bmi_cat_sr_d := cut(bmi_sr_inv, bmi_breaks, bmi_levels, right = FALSE)]
  dt[]
}

chk <- derive_bmi(copy(co)[!is.na(bmi_adj)])
agree_adj <- chk[, mean(bmi_cat_d == bmi_cat)]
agree_sr  <- chk[, mean(bmi_cat_sr_d == bmi_cat_sr)]
cat(sprintf("Category agreement on observed data: adjusted %.5f, self-reported %.5f\n", agree_adj, agree_sr))
stopifnot(agree_adj == 1, agree_sr > 0.999)

# Imputation data -------------------------------------------------------------

imp_vars <- c("htn", "dm", "hd", "bmi_adj", "educ", "income", "smoking", "alcohol", "phys_act")
dat <- co[, .(ADM_RNO, sex, age_grp, age_adj, province, logw = log(WTS_M),
              htn = factor(htn, 0:1), dm = factor(dm, 0:1), hd = factor(hd, 0:1),
              bmi_adj,
              educ = factor(educ, ordered = TRUE),
              income = factor(income, ordered = TRUE),
              phys_act = factor(phys_act, ordered = TRUE),
              smoking, alcohol)]
dat[, stratum := interaction(sex, age_grp, sep = " ", drop = TRUE)]

meth <- make.method(as.data.frame(dat[, !c("ADM_RNO", "sex", "age_grp", "stratum")]))
meth[c("htn", "dm", "hd")] <- "logreg"
meth["bmi_adj"] <- "pmm"
meth[c("educ", "income", "phys_act")] <- "polr"
meth[c("smoking", "alcohol")] <- "polyreg"
stopifnot(all(meth[setdiff(names(meth), imp_vars)] == ""))

run_one <- function(s, i) {
  d <- droplevels(as.data.frame(dat[stratum == s, !c("ADM_RNO", "sex", "age_grp", "stratum")]))
  mt <- meth[names(d)]
  pred <- make.predictorMatrix(d)
  mice(d, m = 1, maxit = maxit, method = mt, predictorMatrix = pred,
       seed = seed + 1000L * match(s, strata) + i, printFlag = FALSE)
}

strata <- levels(dat$stratum)
jobs <- CJ(s = strata, i = seq_len(m), sorted = FALSE)
t0 <- Sys.time()
fits <- mclapply(seq_len(nrow(jobs)), function(j) run_one(jobs$s[j], jobs$i[j]),
                 mc.cores = n_core, mc.preschedule = FALSE)
cat("Imputation time:", format(Sys.time() - t0), "\n")
err <- vapply(fits, inherits, logical(1), "try-error")
if (any(err)) stop("Imputation failed: ", paste(unique(unlist(fits[err])), collapse = "; "))

mids_by_stratum <- lapply(strata, function(s) {
  idx <- which(jobs$s == s)
  Reduce(ibind, fits[idx])
})
names(mids_by_stratum) <- strata
logged <- rbindlist(lapply(strata, function(s) {
  lg <- mids_by_stratum[[s]]$loggedEvents
  if (is.null(lg)) NULL else data.table(stratum = s, lg)
}), fill = TRUE)
saveRDS(mids_by_stratum, file.path(derived_dir, "mids_by_stratum.rds"))

# Completed data sets ---------------------------------------------------------

to_int <- function(f) as.integer(as.character(f))

imputed <- lapply(seq_len(m), function(k) {
  dk <- rbindlist(lapply(strata, function(s) {
    cmp <- as.data.table(complete(mids_by_stratum[[s]], k))
    cmp[, ADM_RNO := dat[stratum == s]$ADM_RNO]
    cmp[, .(ADM_RNO, htn, dm, hd, bmi_adj, educ, income, smoking, alcohol, phys_act)]
  }))
  dk <- merge(co[, .(ADM_RNO, WTS_M, sex, age_grp, age_adj, province, proxy, bmi_cat_obs = bmi_cat)],
              dk, by = "ADM_RNO", sort = FALSE)
  dk[, `:=`(htn = to_int(htn), dm = to_int(dm), hd = to_int(hd),
            educ = factor(educ, ordered = FALSE), income = factor(income, ordered = FALSE),
            phys_act = factor(phys_act, ordered = FALSE))]
  derive_bmi(dk)
  dk[, `:=`(bmi_cat = bmi_cat_d, bmi_cat_sr = bmi_cat_sr_d, bmi_sr = bmi_sr_inv,
            bmi_cat_d = NULL, bmi_cat_sr_d = NULL, bmi_sr_inv = NULL)]
  dk[, n_cond := htn + dm + hd]
  dk[, cmm := as.integer(n_cond >= 2L)]
  dk[, bmi_imputed := is.na(bmi_cat_obs)][, bmi_cat_obs := NULL]
  dk[, .imp := k]
  setorder(dk, ADM_RNO)
  stopifnot(!anyNA(dk[, .(htn, dm, hd, bmi_adj, educ, income, smoking, alcohol, phys_act, cmm)]))
  dk[]
})
saveRDS(imputed, file.path(derived_dir, "imputed.rds"))

# Diagnostics -----------------------------------------------------------------

# Trace of chain means by iteration, per stratum.
trace <- rbindlist(lapply(strata, function(s) {
  cm <- mids_by_stratum[[s]]$chainMean
  rbindlist(lapply(dimnames(cm)[[1]], function(v) {
    x <- cm[v, , , drop = TRUE]
    if (all(is.na(x))) return(NULL)
    data.table(stratum = s, variable = v, iteration = rep(seq_len(dim(cm)[2]), dim(cm)[3]),
               chain = rep(seq_len(dim(cm)[3]), each = dim(cm)[2]), mean = as.vector(x))
  }))
}))
p_trace <- ggplot(trace, aes(iteration, mean, group = chain)) +
  geom_line(alpha = 0.35, linewidth = 0.3) +
  facet_grid(variable ~ stratum, scales = "free_y") +
  labs(x = "Iteration", y = "Mean of imputed values (factor codes for categorical)") +
  theme_bw(base_size = 7)
ggsave(file.path(fig_dir, "mi_trace.png"), p_trace, width = 11, height = 9, dpi = 200, bg = "white")

# Observed versus imputed distributions (unweighted), pooled over imputations.
all_imp <- rbindlist(imputed)
obs_flags <- co[, .(ADM_RNO,
                    htn_m = is.na(htn), dm_m = is.na(dm), hd_m = is.na(hd), bmi_cat_m = is.na(bmi_cat),
                    educ_m = is.na(educ), income_m = is.na(income), smoking_m = is.na(smoking),
                    alcohol_m = is.na(alcohol), phys_act_m = is.na(phys_act))]
all_imp <- merge(all_imp, obs_flags, by = "ADM_RNO")
cmp_vars <- c("htn", "dm", "hd", "bmi_cat", "educ", "income", "smoking", "alcohol", "phys_act")
ovi <- rbindlist(lapply(cmp_vars, function(v) {
  x <- all_imp[, .(sex, value = as.character(get(v)), imputed = get(paste0(v, "_m")))]
  x[, .N, by = .(sex, imputed, value)][, pct := round(100 * N / sum(N), 1), by = .(sex, imputed)][, variable := v]
}))
ovi[, source := fifelse(imputed, "Imputed", "Observed")]
ovi <- dcast(ovi, variable + sex + value ~ source, value.var = "pct", fill = 0)
fwrite(ovi, file.path(tab_dir, "mi_observed_vs_imputed.csv"))

# Pooled prevalence of multimorbidity by sex and BMI (Rubin's rules).
bsw <- readRDS(file.path(derived_dir, "bsw.rds"))
rw  <- rep_weights(imputed[[1]], bsw)
rm(bsw); invisible(gc())
stopifnot(all(vapply(imputed, function(x) identical(x$ADM_RNO, imputed[[1]]$ADM_RNO), logical(1))))
per_imp <- rbindlist(lapply(imputed, function(dk) {
  rbind(prevalence(dk, rw, "cmm", "sex"), prevalence(dk, rw, "cmm", c("sex", "bmi_cat")), fill = TRUE)[, .imp := dk$.imp[1]]
}))
per_imp[is.na(bmi_cat), bmi_cat := "All"]
pooled <- per_imp[, {
  qbar <- mean(p); u <- mean(se^2); b <- var(p); tv <- u + (1 + 1 / m) * b
  r  <- (1 + 1 / m) * b / u
  df <- (m - 1) * (1 + 1 / r)^2
  se_t <- sqrt(tv)
  .(p = qbar, se = se_t, fmi = (r + 2 / (df + 3)) / (r + 1),
    lcl = plogis(qlogis(qbar) - qt(0.975, df) * se_t / (qbar * (1 - qbar))),
    ucl = plogis(qlogis(qbar) + qt(0.975, df) * se_t / (qbar * (1 - qbar))),
    n_cases = mean(n_num), n = mean(n_den))
}, by = .(sex, bmi_cat)]
pooled[, cv := se / p]
pooled[, flag := release_flag(pmin(n_cases, n), cv)]
pooled[, `:=`(sex = factor(sex, levels(co$sex)), bmi_cat = factor(bmi_cat, c("All", bmi_levels)))]
setorder(pooled, sex, bmi_cat)
fwrite(pooled[, .(sex, bmi_cat, n_cases = round(n_cases), n,
                  prevalence_pct = round(100 * p, 1), lcl_pct = round(100 * lcl, 1),
                  ucl_pct = round(100 * ucl, 1), cv = round(cv, 3), fmi = round(fmi, 3), flag)],
       file.path(tab_dir, "prevalence_cmm_mi.csv"))

print(pooled[, .(sex, bmi_cat, prev = fmt_prev(p, lcl, ucl, flag), fmi = round(fmi, 3))])
cat("Logged events:", if (nrow(logged)) nrow(logged) else 0L, "\n")
if (nrow(logged)) print(unique(logged[, .(stratum, meth, out)]))

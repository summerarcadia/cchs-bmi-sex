# Methods notes

Working notes for the Methods section. Each item names the script that implements it.

## Data and population (`R/01_import.R`, `R/02_cohort.R`)

- CCHS 2017-2018 Annual Component PUMF, 113,290 respondents aged 12 and older.
- Eligible cohort: adults 18 and older, excluding the 25 who were currently pregnant (BMI not derived for pregnant women). n = 104,611.
- Variable definitions: `docs/data_dictionary.md`.

## Survey design and variance (`R/functions.R`)

- Estimates are weighted with the master weight `WTS_M`.
- Variance uses the 1,000 PUMF bootstrap replicate weights (`BSW1`-`BSW1000`), as recommended in the CCHS 2017-2018 User Guide (section 10.3) and PUMF Quality Assurance Guidelines. The PUMF replicate weights are produced with a method that differs from the master and share files for confidentiality (PUMF Complement User Guide, section 3).
- Variance estimator: v = (1/B) Σ (θ_b − θ̂)², B = 1,000, centred on the full-sample estimate. For model fitting: `survey::svrepdesign(type = "bootstrap", combined.weights = TRUE, scale = 1/1000, rscales = 1, mse = TRUE)`.
- Domain estimates for descriptive tables are computed directly as ratios of weighted sums over the full-sample and replicate weights. This is algebraically the same estimator as `survey::svymean` on the replicate design; the two agree to at least 1e-7 on test estimates.
- 95% CIs for proportions: Wald interval on the logit scale (delta method).

## Release rules (`R/functions.R`)

From the CCHS 2017-2018 PUMF Quality Assurance Guidelines:

- Not released (S): estimate based on fewer than 30 respondents. For prevalence, the smaller of the numerator and the denominator is used.
- E: CV 15.1-35.0%, released with a caution.
- F: CV above 35.0%, not released.
- Proportions rounded to one decimal place; population totals rounded to the nearest thousand.

## Descriptive analysis (`R/04_descriptives.R`)

- Table 1: unweighted n and weighted percentage of each characteristic by sex, overall and by adjusted BMI category, on observed data with missing values shown as a category.
- Table 2: weighted prevalence (95% CI) of hypertension, diabetes, heart disease and multimorbidity overall, by sex, by sex and BMI category, and by sex and age group. By sex, age group and BMI category together: `output/tables/prevalence.csv`.

## Missing data (`R/05_impute.R`)

- Missing values in the eligible cohort (unweighted): BMI 6.3% overall (men 4.6%, women 7.8%); education 3.1%; physical activity 2.8%; the conditions, income, smoking and alcohol each under 0.5%. Proxy interviews (n = 2,644) were not asked height and weight. 91,978 adults (87.9%) have complete data on all analysis variables.
- Multiple imputation by chained equations (`mice` 3.19.0), 20 imputations, 10 iterations.
- Imputed separately within each of six sex by age group strata (18-44, 45-64, 65+), so that associations that vary by sex or age, including between BMI and the conditions, are preserved.
- Imputed variables and methods: hypertension, diabetes, heart disease (logistic regression); continuous adjusted BMI (predictive mean matching); education, household income, physical activity (proportional odds); smoking, alcohol (multinomial logistic). Complete predictors: age group (6 levels), province, log survey weight.
- Derived after imputation: adjusted BMI category from adjusted BMI using the Health Canada cut-points; self-reported BMI by inverting the sex-specific correction equation, then its category (agreement with `HWTDGISW` on observed data: 99.987%); multimorbidity from the three imputed conditions.
- Pooling: Rubin's rules, with the bootstrap variance as the within-imputation variance.
- Diagnostics: trace plots of chain means (`output/figures/mi_trace.png`); observed against imputed distributions (`output/tables/mi_observed_vs_imputed.csv`); pooled multimorbidity prevalence with fraction of missing information (`output/tables/prevalence_cmm_mi.csv`).
- Diagnostic results (m = 20, maxit = 10): chains mix with no trend across iterations in all strata; no logged events. Imputed values show higher condition prevalence than observed values (heart disease 15.2% vs 7.9% in men, unweighted), consistent with who is missing: 61.9% of those missing heart disease are 65 or older, against 30.8% of those with it observed. Condition data are missing for fewer than 450 adults each. Fraction of missing information for pooled multimorbidity prevalence by sex and BMI category: 0.01-0.20.

## Regression models (`R/06_models.R`, `R/07_model_tables.R`)

- Survey-weighted logistic regression. Primary outcome: cardiometabolic multimorbidity. Secondary outcomes: hypertension, diabetes and heart disease, each modelled separately.
- Exposure: adjusted BMI category (reference: normal weight). Effect modifier: sex (reference: male).
- Adjustment set, as specified in the protocol: age group (6 levels), household income (5), household education (3), smoking (3), alcohol (3), physical activity (4).
- Two models per outcome: main effects (BMI + sex + covariates), and the same with a sex by BMI category product term.
- Reported from the interaction model: OR for each BMI category against normal weight within men and within women; ratio of ORs (women / men) as the multiplicative interaction; OR for women against men within each BMI category; joint ORs for each sex and BMI combination against normal-weight men (`models_interaction.csv`).
- Weights were rescaled to mean 1 before fitting, as `survey::svyglm` does; with population-scale weights the iterations otherwise fail to converge.
- Fitting: Newton-Raphson on data collapsed to unique covariate pattern and outcome combinations, with full-sample and replicate weights summed within each pattern. This leaves the weighted likelihood unchanged. Each replicate fit starts from the full-sample estimates. Checked against `survey::svyglm` on the replicate design for the multimorbidity interaction model in one imputation: coefficients agree to 1.5e-8 and standard errors to within a ratio of 1.000002.
- Variance: bootstrap covariance of replicate coefficients within each imputation, (1/B) Σ (b_r − b)(b_r − b)', pooled across 20 imputations with Rubin's rules and Rubin's degrees of freedom.
- Global test of the sex by BMI interaction: pooled multivariate Wald test of the three product terms (D1; Li, Raghunathan and Rubin 1991).
- Multiple testing: Holm adjustment across the four global interaction tests (one per outcome); Benjamini-Hochberg q-values across the twelve category-specific ratio-of-OR tests (three BMI categories by four outcomes).
- Estimates for a BMI category resting on fewer than 30 unweighted cases in either sex are flagged and not interpreted (applies to underweight).

## Points for later sessions

- Multimorbidity among underweight adults rests on 13 men and 20 women (Table 2, not released). Odds ratios and additive interaction measures for underweight will be unstable; the additive interaction analysis should focus on overweight and obese against normal weight, with underweight reported in the models but not interpreted.

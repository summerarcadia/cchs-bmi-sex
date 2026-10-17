# Sex differences in BMI and cardiometabolic multimorbidity in Canadian adults

Cross-sectional analysis of the Canadian Community Health Survey (CCHS) 2017-2018 Public Use Microdata File. The study estimates whether the association between BMI category and cardiometabolic multimorbidity (two or more of hypertension, diabetes and heart disease) differs between women and men, on both the multiplicative and additive scales.

## Data

The microdata are not included in this repository. Download them from Statistics Canada:

- Product: Canadian Community Health Survey: Public Use Microdata File, 2000 to 2022 (catalogue 82M0013X2024001)
- Page: https://www150.statcan.gc.ca/n1/pub/82m0013x/82m0013x2024001-eng.htm
- File: `2017-2018.zip`

Unzip it into `data/raw/2017-2018/`, then unzip `Bootstrap weights_Poids bootstrap.zip` into `data/raw/2017-2018/3226_CCHS-Ann_ESCC_2017-2018/bsw/` and its `data.zip` in place. The expected paths are set in `R/00_setup.R`.

## Running

From the project root:

```
Rscript R/01_import.R        # read PUMF variables and bootstrap weights
Rscript R/02_cohort.R        # cohort, recodes, outcome, flow counts
Rscript R/03_flow_diagram.R  # participant flow figure
Rscript R/04_descriptives.R  # Table 1 and weighted prevalence
Rscript R/05_impute.R        # multiple imputation (about 40 minutes on 8 cores)
Rscript R/06_models.R        # survey-weighted logistic models (about 25 minutes on 8 cores)
Rscript R/07_model_tables.R  # odds ratio tables and interaction tests
Rscript R/08_additive.R      # RERI, AP, SI and standardized prevalence differences
Rscript R/09_age_strata.R    # interaction by age group (about 15 minutes)
Rscript R/10_sensitivity.R   # complete case, BMI definitions, underweight excluded (about 30 minutes)
Rscript R/11_propensity.R    # propensity score weighting and balance plots
Rscript R/12_evalues.R       # E-values
Rscript R/13_additive_tables.R  # formatted Session 4 tables
```

`R/functions.R` holds the survey design and release-rule helpers. Methods details are in `docs/methods_notes.md`.

R 4.4. Packages: data.table, readr, here, ggplot2, survey, srvyr, mice, mitools, gtsummary, WeightIt, cobalt, EValue.

## Layout

```
R/         analysis scripts, run in numeric order
docs/      data dictionary and feasibility notes
output/    tables and figures
data/      raw and derived microdata (not tracked)
```

## Data licence and attribution

Adapted from Statistics Canada, Canadian Community Health Survey, 2017-2018 Public Use Microdata File. This does not constitute an endorsement by Statistics Canada of this product. Used under the Statistics Canada Open Licence.

# Sex differences in BMI and cardiometabolic multimorbidity in Canadian adults

Cross-sectional analysis of the Canadian Community Health Survey (CCHS) 2017-2018 Public Use Microdata File. The study estimates whether the association between BMI category and cardiometabolic multimorbidity (two or more of hypertension, diabetes and heart disease) differs between women and men, on both the multiplicative scale (ratio of odds ratios) and the additive scale (RERI and model-standardized prevalence differences), overall and by age group.

The manuscript is built from `manuscript/manuscript.md` into `output/report/manuscript.docx`; every number in it is filled in from the result tables in `output/tables`.

## Data

The microdata are not included in this repository. Download them from Statistics Canada:

- Product: Canadian Community Health Survey: Public Use Microdata File, 2000 to 2022 (catalogue 82M0013X2024001)
- Page: https://www150.statcan.gc.ca/n1/pub/82m0013x/82m0013x2024001-eng.htm
- File: `2017-2018.zip`

Unzip it into `data/raw/2017-2018/`. Inside `3226_CCHS-Ann_ESCC_2017-2018/`, unzip `Bootstrap weights_Poids bootstrap.zip` into a folder named `bsw/`, then unzip its `data.zip` in place. The expected paths are set in `R/00_setup.R`. Some archive members have Latin-1 file names; on macOS, `ditto -x -k` extracts them where `unzip` fails.

## Running

From the project root, in order:

```
Rscript R/01_import.R           # read PUMF variables and bootstrap weights
Rscript R/02_cohort.R           # cohort, recodes, outcome, flow counts
Rscript R/03_flow_diagram.R     # participant flow figure
Rscript R/04_descriptives.R     # Table 1 and weighted prevalence
Rscript R/05_impute.R           # multiple imputation (about 40 minutes on 8 cores)
Rscript R/06_models.R           # survey-weighted logistic models (about 55 minutes)
Rscript R/07_model_tables.R     # odds ratio tables and interaction tests
Rscript R/08_additive.R         # RERI, AP, SI and standardized prevalence differences
Rscript R/09_age_strata.R       # interaction by age group (about 15 minutes)
Rscript R/10_sensitivity.R      # complete case, BMI definitions, underweight excluded (about 20 minutes)
Rscript R/11_propensity.R       # propensity score weighting and balance plots (about 15 minutes)
Rscript R/12_evalues.R          # E-values
Rscript R/13_additive_tables.R  # formatted additive, age, sensitivity and E-value tables
Rscript R/14_figures.R          # report figures
Rscript R/15_report.R           # manuscript and one-page summary (.docx)
```

Scripts 09-11 save their model fits; set `REUSE_FITS=1` to rebuild their tables without refitting. `R/functions.R` holds the survey design, model fitting, pooling and additive interaction helpers. Methods details are in `docs/methods_notes.md`.

Requirements: R 4.4 with data.table, readr, here, ggplot2, survey, mice, gt, WeightIt, cobalt, EValue, officer and flextable.

## Layout

```
R/            analysis scripts, run in numeric order
manuscript/   manuscript and summary text with placeholders for results
docs/         data dictionary, methods notes, feasibility notes, STROBE checklist, statements
output/       result tables, figures and the built documents
data/         raw and derived microdata (not tracked)
```

## Licence

Code: MIT licence (see `LICENSE`).

Data: adapted from Statistics Canada, Canadian Community Health Survey, 2017-2018 Public Use Microdata File. This does not constitute an endorsement by Statistics Canada of this product. Used under the Statistics Canada Open Licence.

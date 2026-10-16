# Session 1 feasibility check

Date: 2026-10-13
Source files: Statistics Canada, CCHS Public Use Microdata File, 2000 to 2022 (82M0013X2024001), Statistics Canada Open Licence.
Cycles checked: 2022 (67,079 records) and 2019-2020 (108,252 records). Variable names, universes and codes below were read from each cycle's `DataDictionary_Freqs` PDF, the 2022 PUMF Complement User Guide and the 2022 PUMF Grouped Variables document, then confirmed against the CSV files.

## Status: stopped, key variables missing

Two variables needed by the protocol are not available in the form the protocol requires, in either cycle.

### 1. BMI category (exposure)

| Cycle | Variable | Categories on PUMF |
|---|---|---|
| 2022 | `HWTDGBCC` (adjusted), `HWTDGISW` (self-reported) | 1 = underweight/normal, 2 = overweight/obese I-III |
| 2019-2020 | `HWTDGBCC`, `HWTDGISW` | same two categories |

The full classification (`HWTDVBCC`, `HWTDGCOR`) exists on the master file only. The PUMF Complement User Guide lists BMI status among the variables regrouped for disclosure control.

Consequences: no separate obese category, no underweight category, and none of the planned sensitivity analyses on cut-points or excluding underweight can be run on the PUMF.

In 2022, 452 adults have `HWTDGBCC = 6` (valid skip). The derivation excludes pregnant women, so these are most likely pregnant respondents; pregnancy itself is suppressed on the PUMF. To confirm before treating them as an exclusion.

### 2. Heart disease (outcome component)

There is no heart disease item on its own in either PUMF. The only cardiovascular variable is `CCCDGCAR`, which groups heart disease, past heart attack and the effects of stroke.

| Cycle | Universe |
|---|---|
| 2022 | Respondents aged 35 and older (ages 18-34 are valid skip) |
| 2019-2020 | All respondents |

Proposed substitute: `CCCDGCAR`, labelled "heart disease or stroke" throughout.

### 3. Physical activity (covariate)

| Cycle | Inclusion flag | Adults with the module |
|---|---|---|
| 2022 | `DOPAA` | about 1,800 of 67,079 (`PAADVWHO`, `PAADVAC2`) |
| 2019-2020 | `DOPAA` | 22,006 of 108,252 |

Physical activity is optional content, selected by only some provinces, so it is missing by design for most respondents. Multiple imputation is not appropriate for data missing by design at this scale.

## Variables that are available

| Concept | 2022 | 2019-2020 | Notes |
|---|---|---|---|
| Sex | `DHH_SEX` | `DHH_SEX` | 2022 label: sex at birth |
| Age group | `DHHGAGE` | `DHHGAGE` | 12-17, 18-34, 35-49, 50-64, 65+ |
| Hypertension | `CCC_80` | `CCC_065` | "Do you have high blood pressure?" |
| Diabetes | `CCC_05` | `CCC_095` | "Do you have diabetes?" |
| Heart disease or stroke | `CCCDGCAR` | `CCCDGCAR` | see above |
| Household income | `INCDGHH`, `INCDGRCA` | `INCDGHH`, `INCDGRCA` | 5 groups; national quintile |
| Education (household) | `EDDVH3` | `EHG2DVH3` | 3 levels; renamed between cycles |
| Smoking | `SMKDVSTY` | `SMKDVSTY` | 6 categories |
| Alcohol | `ALCDVTTM` | `ALCDVTTM` | regular / occasional / none in 12 months |
| Province | `GEOGPRV` | `GEOGPRV` | |
| Survey weight | `WTS_M` | `WTS_M` | |
| Bootstrap weights | `BSW1`-`BSW1000` in `cchs_escc_bsw.csv` | same | linked on `ADM_RNO` (2022) / `ADM_RNO1` (2019-2020) |

## Variance estimation

Both PUMFs ship with 1,000 bootstrap replicate weights. The 2022 Complement User Guide states they are produced with a different method from the master and share file weights for confidentiality, that all are strictly positive, and that removing one health region (Nord-du-Québec) leaves replicate sums slightly unequal with negligible effect on variance. Planned approach: `survey::svrepdesign(type = "bootstrap")` with the PUMF replicate weights.

## Age strata

The PUMF age groups are 18-34, 35-49, 50-64 and 65+. The protocol's 18-44 / 45-64 split cannot be built.

## Update: 2017-2018 PUMF resolves the gaps

Checked 2026-10-15. File: `CCHS-82M0021X-E-2017-2018-Annual_F1.csv` (113,290 records, 1,051 variables), with the 2017-2018 PUMF Data Dictionary and Complement User Guide.

| Gap in 2019-2020 / 2022 | 2017-2018 PUMF |
|---|---|
| BMI in two groups only | `HWTDGBCC` (adjusted) and `HWTDGISW` (self-reported): underweight, normal, overweight, obese I-III. Continuous BMI also available: `HWTDGCOR` (adjusted, 13.92-57.48) and `HWTDGBMI` (self-reported, 12.07-54.42), top- and bottom-coded |
| No separate heart disease item | `CCC_085` "Do you have heart disease?", all respondents with `DOCCC = 1`; stroke separately as `CCC_090` |
| CVD only for 35+ (2022) | `CCC_085` asked at all ages |
| Physical activity missing by design | `DOPAA = 1` for all 113,290 records; `PAADVAC2` defined for all adults |
| Age only in 4 adult groups | `DHHGAGE` in 5-year groups (18-19, 20-24, ... 80+), so 18-44 / 45-64 / 65+ can be built |

Hypertension `CCC_065`, diabetes `CCC_095`, education `EHG2DVH3`, income `INCDGHH`, smoking `SMKDVSTY`, alcohol `ALCDVTTM`, sex `DHH_SEX` and weight `WTS_M` are all present. The Complement User Guide states that the 2017-2018 PUMF is the first to include bootstrap weights, produced by the same PUMF-specific method described for 2022.

Decision (2026-10-16): 2017-2018 is the analysis cycle. Physical activity is retained as a covariate. The trade-off is data collected before COVID-19, to be stated as a limitation.

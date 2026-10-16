# Data dictionary

Source: Statistics Canada, Canadian Community Health Survey, Annual Component, 2017-2018 Public Use Microdata File (82M0013X2024001), Statistics Canada Open Licence.
Codes are from the 2017-2018 PUMF Data Dictionary and Derived Variables Specifications. Derived analysis variables are built in `R/02_cohort.R`.

Codes 6 (valid skip), 7 (don't know), 8 (refusal) and 9 (not stated), and 999.96 / 999.99 for continuous BMI, are set to missing unless noted.

## Design

| Analysis variable | Source | Definition |
|---|---|---|
| `ADM_RNO` | `ADM_RNO` | Record number; links the PUMF and bootstrap files |
| `WTS_M` | `WTS_M` | Survey weight (equal to `FWGT` in the bootstrap file) |
| Replicate weights | `BSW1`-`BSW1000` | Bootstrap file `bsw.txt`, fixed width, two implied decimals |
| `province` | `GEO_PRV` | Province or territory (13 levels) |
| `proxy` | `ADM_PRX` | Health component completed by proxy (height and weight not asked) |

## Population

| Criterion | Source | Rule |
|---|---|---|
| Adult | `DHHGAGE` | Codes 03-16 (18 and older) |
| Currently pregnant (excluded) | `HWTDGBCC` | Code 6 in an adult. The derivation sets 6 when `MAC_025 = 1`; `MAC_025` is not on the PUMF |

## Exposure

| Analysis variable | Source | Coding |
|---|---|---|
| `bmi_cat` (primary) | `HWTDGBCC` | BMI adjusted for self-report bias: 1 Underweight (<18.5), 2 Normal (18.5-24.99), 3 Overweight (25-29.99), 4 Obese (>=30). Reference: Normal |
| `bmi_cat_sr` | `HWTDGISW` | Same categories, unadjusted self-reported BMI |
| `bmi_adj` | `HWTDGCOR` | Continuous adjusted BMI, grouped by Statistics Canada to 13.92-57.48 |
| `bmi_sr` | `HWTDGBMI` | Continuous self-reported BMI, grouped to 12.07-54.42 |

Adjusted BMI applies sex-specific correction equations to self-reported BMI (Connor Gorber et al. 2008; Shields et al. 2011): males −1.07575 + 1.07592 × BMI, females −0.12374 + 1.05129 × BMI. Because the correction differs by sex, results on adjusted and self-reported BMI will be compared.

## Effect modifier

| Analysis variable | Source | Coding |
|---|---|---|
| `sex` | `DHH_SEX` | 1 Male, 2 Female. Question: "Is [respondent name] male or female?" |

## Outcome and components

| Analysis variable | Source | Coding |
|---|---|---|
| `htn` | `CCC_065` | "Do you have high blood pressure?" 1 = yes, 0 = no |
| `dm` | `CCC_095`, `CCC_105`, `CCC_110` | "Do you have diabetes?" 1 = yes, 0 = no. Recoded to 0 when diabetes was first diagnosed during pregnancy (`CCC_105 = 1`) and never otherwise (`CCC_110 = 2`); missing when `CCC_110` is don't know or not stated |
| `hd` | `CCC_085` | "Do you have heart disease?" 1 = yes, 0 = no |
| `n_cond` | derived | Number of the three conditions reported (0-3) |
| `cmm` | derived | Cardiometabolic multimorbidity: 1 if `n_cond >= 2`; 0 if the known conditions plus missing ones cannot reach 2; otherwise missing |

Known issue (Data Dictionary note on `CCC_105`): in 2018, some females aged 15+ who answered "don't know" or "refused" to age at diabetes diagnosis were coded "not applicable" instead of "not stated" on `CCC_105`. These respondents keep `dm = 1`.

## Covariates

| Analysis variable | Source | Coding |
|---|---|---|
| `age_grp` | `DHHGAGE` | Stratification: 18-44 (03-08), 45-64 (09-12), 65+ (13-16) |
| `age_adj` | `DHHGAGE` | Adjustment: 18-34, 35-44, 45-54, 55-64, 65-74, 75+ |
| `income` | `INCDGHH` | Total household income: <20k, 20-39k, 40-59k, 60-79k, 80k+ |
| `educ` | `EHG2DVH3` | Highest household education: less than secondary, secondary, post-secondary |
| `smoking` | `SMKDVSTY` | Never (05 experimental, 06 lifetime abstainer), Former (03, 04), Current (01, 02) |
| `alcohol` | `ALCDVTTM` | None in past 12 months (3), Occasional (2), Regular (1) |
| `phys_act` | `PAADVWHO` | Weekly activity volume, MET-min: 0 (level 4), 1-449 (3), 450-899 (2), >=900 (1) |

## Not available on the PUMF

| Concept | Note |
|---|---|
| Pregnancy status (`MAC_025`) | Suppressed; inferred from `HWTDGBCC = 6` |
| Age in single years | Grouped in 5-year bands; menopausal status cannot be derived |
| Measured height, weight, blood pressure, glucose | Not collected in CCHS |

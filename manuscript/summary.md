# Obesity and cardiometabolic multimorbidity in Canadian adults: relative and absolute sex differences

Hayley So | h2so@uwaterloo.ca

## The question

Does the link between body weight and having two or more of hypertension, diabetes and heart disease differ between women and men? Most studies compare odds ratios (the relative scale). This study also measures the absolute excess, the extra percentage of people with multimorbidity, which is what matters for planning prevention and services.

## Data and methods

- Canadian Community Health Survey 2017-2018 public use file: {{n_cohort}} adults, weighted to represent the Canadian household population.
- BMI from self-reported height and weight, adjusted for reporting bias. Survey-weighted logistic regression adjusted for age, income, education, smoking, alcohol and physical activity, with bootstrap variance and {{m}} multiple imputations.
- Interaction measured on both scales: ratio of odds ratios, relative excess risk due to interaction (RERI), and model-standardized prevalence differences.

## Main findings

- Obesity was associated with higher odds of multimorbidity in women (odds ratio {{or_cmm_women_obese}}) and men ({{or_cmm_men_obese}}). On the relative scale there was no clear sex difference.
- On the absolute scale, obesity was associated with {{rd_female_cmm_obese_pt}} more percentage points of multimorbidity in women and {{rd_male_cmm_obese_pt}} in men (difference {{rd_diff_cmm_obese}}; RERI {{reri_cmm_obese}}).
- The absolute excess rose steeply with age. At 65 and older it was {{rd_female_obese_65}} percentage points in women and {{rd_male_obese_65}} in men.
- For diabetes alone, the relative association with obesity was stronger in women (ratio of odds ratios {{or_dm_ratio_obese}}), yet the absolute excess was similar in both sexes.
- Findings were consistent with complete cases, unadjusted self-reported BMI and covariate-adjusted propensity score weighting.

## Why it matters

Looking only at relative measures suggests obesity matters as much or more for women. The absolute scale shows more of the obesity-associated burden of multimorbidity among men, especially older men. Both views are needed to decide where prevention will reach the most people.

## Limitations

Cross-sectional data show associations, not cause and effect. Weight, height and conditions were self-reported, so undiagnosed conditions were missed and misclassification is likely. Menopausal status and gender were not available. Data were collected before the COVID-19 pandemic.

## Data and code

Statistics Canada public use microdata (Open Licence). Analysis code: [repository URL to be added]. Preprint: [link to be added].

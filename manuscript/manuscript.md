# Sex differences in the association between body mass index and cardiometabolic multimorbidity in Canadian adults: additive and multiplicative interaction in the Canadian Community Health Survey

Hayley So

[Affiliation: to be confirmed by author]

Correspondence: h2so@uwaterloo.ca

## Abstract

Background: Sex differences in the association between body mass index (BMI) and cardiometabolic disease are usually reported as ratios of odds ratios. Interaction on the additive scale, which reflects absolute excess burden and is more relevant to public health targeting, has not been reported for cardiometabolic multimorbidity in Canadian adults. We estimated sex differences in the association between BMI category and cardiometabolic multimorbidity on both scales.

Methods: Cross-sectional analysis of the Canadian Community Health Survey 2017-2018 Public Use Microdata File ({{n_cohort}} adults). Multimorbidity was two or more of self-reported hypertension, diabetes and heart disease. BMI was self-reported and adjusted for reporting bias. We fitted survey-weighted logistic regression with a sex by BMI interaction, adjusted for age, income, education, smoking, alcohol and physical activity, using bootstrap replicate weights and {{m}} multiple imputations. Additive interaction was estimated as the relative excess risk due to interaction (RERI) and as model-standardized prevalence differences.

Results: Obesity was associated with higher odds of multimorbidity in women (odds ratio {{or_cmm_women_obese}}) and men ({{or_cmm_men_obese}}), with no evidence of multiplicative interaction (ratio of odds ratios {{or_cmm_ratio_obese}}). On the additive scale, obesity was associated with a larger absolute excess in men: the standardized prevalence difference was {{rd_female_cmm_obese}} percentage points in women and {{rd_male_cmm_obese}} in men (difference {{rd_diff_cmm_obese}}), and the RERI for obesity and male sex was {{reri_cmm_obese}}. The absolute excess increased with age, and the sex difference was most apparent at ages 65 and older. Results were similar with complete cases, self-reported BMI and covariate-adjusted propensity score weighting, and attenuated when the sexes were not compared at a common covariate distribution.

Conclusions: The association between obesity and cardiometabolic multimorbidity was similar or slightly stronger in women on the relative scale but larger in men on the absolute scale. Reporting interaction on both scales changes which group carries the larger obesity-associated burden.

## Introduction

Excess body weight and the cardiometabolic conditions it accompanies are common in Canada, and national surveys show rising prevalence of obesity, hypertension and diabetes over the past two decades [1-3]. Canadian Community Health Survey (CCHS) analyses have described the co-occurrence of diabetes with cardiovascular conditions [4], the prevalence and determinants of multimorbidity [5,6], sex and gender differences in hypertension and in diabetes care [7,8], and lifestyle factors associated with cardiovascular disease [9].

When studies ask whether the association between BMI and cardiometabolic disease differs between women and men, they usually compare odds ratios, that is, they test for interaction on the multiplicative scale. Interaction on the additive scale asks a different question: whether the absolute excess burden associated with a risk factor differs between groups. Because women and men differ in their baseline prevalence of cardiometabolic disease, the two scales can disagree, and the additive scale is the one that bears on where the largest number of excess cases occurs [10,11]. Additive interaction between obesity and sex has been examined for arthritis in the CCHS [12] and for the accumulation of chronic conditions in a United States cohort [13], but, to our knowledge, it has not been reported for cardiometabolic multimorbidity in Canadian adults.

We used the CCHS 2017-2018 Public Use Microdata File to estimate whether the association between BMI category and cardiometabolic multimorbidity differs between women and men, on both the multiplicative and the additive scale, overall and by age group. The analysis is descriptive and reports associations; it does not estimate causal effects.

## Methods

### Data source and population

The CCHS is a cross-sectional survey of the household population aged 12 and older in Canada, conducted by Statistics Canada. We used the 2017-2018 Annual Component Public Use Microdata File (PUMF), which includes {{n_all}} respondents and 1,000 bootstrap replicate weights [22]. This cycle was chosen because later PUMFs (2019-2020 and 2022) group BMI into two categories, do not separate heart disease from stroke, and collected physical activity from a minority of respondents. We included adults aged 18 and older and excluded those currently pregnant, for whom BMI is not derived ({{n_preg}}), giving {{n_cohort}} adults (supplementary flow diagram).

### Measures

The exposure was BMI category based on self-reported height and weight, adjusted by Statistics Canada for reporting bias using sex-specific correction equations [14,15]: underweight (under 18.5 kg/m2), normal weight (18.5 to 24.9), overweight (25.0 to 29.9) and obesity (30.0 or above). Hypertension, diabetes and heart disease were self-reported diagnosed conditions; diabetes first diagnosed during pregnancy and not otherwise was not counted. The primary outcome, cardiometabolic multimorbidity, was two or more of these three conditions; each condition was a secondary outcome. Sex was the respondent's reported sex (male or female). Covariates were age group (18-34, 35-44, 45-54, 55-64, 65-74, 75 and older), household income (5 groups), household education (3 levels), smoking status (never, former, current), alcohol use in the past 12 months (none, occasional, regular) and weekly physical activity volume (0, 1-449, 450-899, 900 or more MET-minutes). Variable definitions are in the supplementary data dictionary.

### Missing data

Of the {{n_cohort}} adults, {{n_complete}} ({{pct_complete}}%) had complete data. BMI was missing for {{miss_bmi_men}}% of men and {{miss_bmi_women}}% of women, partly because height and weight were not asked in proxy interviews (n = {{n_proxy}}). We used multiple imputation by chained equations [16] with {{m}} imputations, run separately within each sex and age group (18-44, 45-64, 65+) so that associations that differ by sex and age were preserved. Continuous adjusted BMI was imputed by predictive mean matching, and BMI categories and multimorbidity were derived from the imputed values. Estimates were combined with Rubin's rules.

### Statistical analysis

All estimates used the survey weights. Variance was estimated with the PUMF bootstrap replicate weights, with each imputation's bootstrap variance used as the within-imputation variance. We fitted logistic regression for each outcome with main effects for BMI category and sex and then with a sex by BMI category product term, adjusted for the covariates listed above. Estimates were checked against the survey package in R [17]. We report odds ratios for each BMI category against normal weight within each sex, and their ratio (women over men) as the measure of multiplicative interaction. The global interaction test was a pooled Wald test of the three product terms [18], with Holm adjustment across the four outcomes; category-specific ratios carry Benjamini-Hochberg q-values.

Additive interaction was estimated in two ways. First, on the odds ratio scale, we calculated the RERI, attributable proportion (AP) and synergy index (SI), using normal-weight women, the lowest-risk group, as the reference so that the second exposure is male sex [19]. Second, because odds ratios approximate risk ratios poorly for common outcomes, we estimated model-standardized prevalence: predicted probabilities for each sex and BMI category were averaged over the survey-weighted covariate distribution of the whole cohort. The BMI-associated prevalence difference was calculated within each sex, and the difference between women and men is the interaction contrast on the probability scale. For both approaches, variance was obtained by repeating the calculation on each of the 1,000 bootstrap replicates within each imputation and pooling across imputations.

We repeated the analysis within age groups (18-44, 45-64, 65+), excluding underweight adults, and conducted sensitivity analyses using complete cases, self-reported (unadjusted) BMI, obesity split into class I (30.0 to 34.9) and class II-III (35.0 or above), excluding underweight adults, and propensity score weighting for BMI category within each sex, with and without covariate adjustment in the outcome model. We calculated E-values for the sex-specific odds ratios [20]. Estimates based on fewer than 30 unweighted cases were not interpreted, following Statistics Canada release guidelines. Underweight results are shown but not interpreted because multimorbidity among underweight adults rested on {{n_uw_cases_m}} cases in men. Reporting follows STROBE [21]. Analyses used R 4.4.2 with the survey, mice, WeightIt, cobalt and EValue packages; code is publicly available.

## Results

The {{n_cohort}} adults represented the Canadian household population aged 18 and older (Table 1). The weighted prevalence of cardiometabolic multimorbidity was {{prev_cmm}}: {{prev_cmm_m}}% in men and {{prev_cmm_f}}% in women. Among adults with normal weight it was {{prev_cmm_m_normal}}% in men and {{prev_cmm_f_normal}}% in women, and among adults with obesity {{prev_cmm_m_obese}}% and {{prev_cmm_f_obese}}%.

[[table:table1]]

Table 1. Characteristics of adults by sex, CCHS 2017-2018. Unweighted n (weighted %). E: coefficient of variation 15.1-35.0%.

### Multiplicative scale

Overweight and obesity were associated with higher odds of multimorbidity in both sexes (Table 2, Figure 1). For obesity, the adjusted odds ratio was {{or_cmm_women_obese}} in women and {{or_cmm_men_obese}} in men; the ratio of odds ratios was {{or_cmm_ratio_obese}} (global interaction {{pint_cmm}}; Holm-adjusted {{pholm_cmm}}). Among the secondary outcomes, the association with diabetes was stronger in women for obesity (ratio {{or_dm_ratio_obese}}) and overweight (ratio {{or_dm_ratio_overweight}}; global interaction Holm-adjusted {{pholm_dm}}). There was no evidence of multiplicative interaction for hypertension (Holm-adjusted {{pholm_htn}}) or heart disease (Holm-adjusted {{pholm_hd}}).

[[table:table2]]

Table 2. Odds ratios (95% CI) for overweight and obesity against normal weight in the same sex, from survey-weighted logistic regression with a sex by BMI interaction, adjusted for age group, household income, household education, smoking, alcohol and physical activity. Interaction p: pooled Wald test of the three product terms (including underweight); Holm adjustment across the four outcomes.

[[figure:figure1_forest.png|6.3|6.5]]

Figure 1. Sex-specific odds ratios for overweight and obesity against normal weight, for each outcome.

### Additive scale

On the additive scale, obesity was associated with a larger absolute excess of multimorbidity in men than in women (Table 3). The standardized prevalence difference for obesity against normal weight was {{rd_female_cmm_obese}} percentage points in women and {{rd_male_cmm_obese}} in men; the difference, women minus men, was {{rd_diff_cmm_obese}}. The RERI for obesity and male sex was {{reri_cmm_obese}}, the AP {{ap_cmm_obese}} and the SI {{si_cmm_obese}}. For overweight, the corresponding difference was {{rd_diff_cmm_overweight}} and the RERI {{reri_cmm_overweight}}. For heart disease, the absolute excess associated with obesity was also larger in men (difference {{rd_diff_hd_obese}}). For diabetes, the larger relative association in women did not translate into a larger absolute excess: the prevalence difference, women minus men, was {{rd_diff_dm_obese}}.

[[table:table3]]

Table 3. Additive interaction between BMI category and sex. Standardized prevalence difference: predicted prevalence from the interaction model averaged over the survey-weighted covariate distribution of the cohort, for each BMI category against normal weight; a negative difference (women minus men) indicates a larger absolute excess in men. RERI and AP use normal-weight women as the reference, so the second exposure is male sex; RERI above 0 indicates a more than additive joint association. Estimates (95% CI) pooled over {{m}} imputations.

### Age groups

The absolute excess of multimorbidity associated with obesity increased with age, from {{rd_female_obese_1844}} percentage points in women and {{rd_male_obese_1844}} in men aged 18-44 to {{rd_female_obese_4564}} and {{rd_male_obese_4564}} at 45-64 and {{rd_female_obese_65}} and {{rd_male_obese_65}} at 65 and older (Figure 2). The difference between women and men was {{rd_diff_obese_4564}} at 45-64 and {{rd_diff_obese_65}} at 65 and older, where the RERI was {{reri_obese_65}}. Estimates for ages 18-44 rest on {{n_ref_cases_young_m}} and {{n_ref_cases_young_f}} cases among normal-weight men and women and are not interpreted.

[[figure:figure2_additive.png|6.3|4.0]]

Figure 2. Standardized difference in multimorbidity prevalence associated with overweight and obesity, by sex and age group.

### Sensitivity analyses

The RERI for obesity was {{reri_complete_case}} with complete cases, {{reri_self_report}} with self-reported BMI and {{reri_excl_underweight}} with underweight adults excluded. Additive interaction was concentrated in obesity class II-III (RERI {{reri_class2}}; class I {{reri_class1}}). Propensity score weighting balanced covariates across BMI categories within each sex (maximum standardized mean difference {{smd_unadj}} before and {{smd_wtd}} after weighting). With covariate adjustment in the weighted outcome model, the RERI was {{reri_ps_dr}} and the prevalence difference {{rd_diff_ps_dr}}. Without covariate adjustment, which compares each sex at its own covariate distribution, the RERI was {{reri_ps_marg}} and the prevalence difference {{rd_diff_ps_marg}}. E-values for the obesity odds ratios for multimorbidity were {{ev_cmm_obese_f}} in women and {{ev_cmm_obese_m}} in men (CI limits {{ev_cmm_obese_f_ci}} and {{ev_cmm_obese_m_ci}}).

## Discussion

In a large national sample of Canadian adults, obesity was strongly associated with cardiometabolic multimorbidity in both sexes. On the relative scale the association was similar, or slightly stronger in women, but on the absolute scale it was larger in men. The additive interaction grew with age and was most evident at 65 and older, and it was concentrated in class II-III obesity.

The disagreement between scales follows from the difference in baseline risk. Normal-weight women had lower adjusted odds of multimorbidity than normal-weight men (odds ratio {{or_female_normal}}). A similar or larger relative association applied to a lower baseline in women and a higher baseline in men produces a larger absolute excess in men. Analyses that report only ratios of odds ratios would conclude that the obesity association is, if anything, stronger in women; the additive scale shows that more of the obesity-associated burden of multimorbidity falls on men. The pattern for diabetes illustrates the same point from the other direction: a clearly stronger relative association in women coexisted with a similar absolute excess in both sexes. The RERI for diabetes, which is computed from odds ratios, pointed towards a larger excess in men, while the standardized prevalence difference did not; for outcomes that are common in some subgroups, measures on the probability scale are the safer basis for statements about absolute excess.

The sex difference on the absolute scale held with complete cases, with unadjusted self-reported BMI and with covariate-adjusted propensity score weighting. It was attenuated when the sexes were compared each at their own covariate distribution without adjustment, which indicates that part of the difference depends on comparing women and men at a common covariate profile, for example the same age distribution. Previous CCHS work found a more than additive joint association of obesity and female sex for arthritis [12], the opposite direction to our finding for cardiometabolic multimorbidity, which suggests that the sex pattern of obesity-associated excess burden differs between chronic conditions.

Strengths of this study include a large, nationally representative sample, design-based variance using the survey's bootstrap replicate weights throughout, multiple imputation that preserves sex and age differences, and the reporting of interaction on both scales with absolute measures that do not rely on odds ratios approximating risk ratios.

There are important limitations. The data are cross-sectional, so associations cannot be interpreted as cause and effect, and the direction of association between BMI and conditions cannot be established. BMI and all three conditions were self-reported, so misclassification is likely. Self-reported height and weight underestimate BMI, and the correction equations used to adjust them differ by sex [14,15]; residual sex-differential misclassification could bias the interaction estimates, although results with unadjusted self-reported BMI were similar. Without measured blood pressure or glucose, undiagnosed hypertension and diabetes were missed, and differences in diagnosis between women and men could contribute to the observed patterns. Age was available only in five-year groups and menopausal status was not available, so menopause could not be examined. Sex was recorded as male or female and gender was not measured. The data were collected before the COVID-19 pandemic. Multimorbidity among underweight adults and among adults aged 18-44 was too rare for reliable estimates. Unmeasured confounding cannot be excluded, although an unmeasured confounder would need associations of at least {{ev_cmm_obese_m_ci}} (men) and {{ev_cmm_obese_f_ci}} (women) on the risk ratio scale with both obesity and multimorbidity to move the confidence intervals of those odds ratios to include 1; E-values were smaller for some secondary associations (supplementary Table S3). In the propensity score analysis, weights were held fixed across bootstrap replicates, so its confidence intervals do not include uncertainty in the propensity model.

In conclusion, the association between obesity and cardiometabolic multimorbidity in Canadian adults was similar on the relative scale in women and men but larger in men on the absolute scale, especially at older ages. Reporting both scales is needed to describe where the obesity-associated burden of multimorbidity is greatest.

## Declarations

Ethics: This study used only the CCHS Public Use Microdata File, which Statistics Canada releases publicly under the Statistics Act and the Statistics Canada Open Licence after removing identifiers and applying disclosure control. Under the Tri-Council Policy Statement: Ethical Conduct for Research Involving Humans (TCPS 2, 2022), research that relies exclusively on publicly available information of this kind (Article 2.2) or on secondary use of anonymous information (Article 2.4) does not require research ethics board review. No ethics review was sought.

Data availability: The CCHS 2017-2018 Public Use Microdata File and its documentation are freely available from Statistics Canada (catalogue 82M0013X2024001, https://www150.statcan.gc.ca/n1/pub/82m0013x/82m0013x2024001-eng.htm) under the Statistics Canada Open Licence. The analysis code is available at https://github.com/summerarcadia/cchs-bmi-sex. Adapted from Statistics Canada, Canadian Community Health Survey, 2017-2018 Public Use Microdata File. This does not constitute an endorsement by Statistics Canada of this product.

Supplementary material: The flow diagram, data dictionary, supplementary Tables S1 (age groups), S2 (sensitivity analyses) and S3 (E-values), covariate balance plots and imputation diagnostics are available in the code repository.

Funding: [to be confirmed by author]

Competing interests: [to be confirmed by author]

## References

1. Dai H, Younis A, Kong JD, Bragazzi NL, Wu J. Trends and regional variation in prevalence of cardiovascular risk factors and association with socioeconomic status in Canada, 2005-2016. JAMA Netw Open. 2021;4(8):e2121443. doi:10.1001/jamanetworkopen.2021.21443

2. Chen KYA, Eddeen AB, Bennett C, Yusuf W, Hennessy D, Barnes JD, et al. Trends in cardiovascular risk factors in Canada: variation by migration and temporal factors, 2001-2018. CJC Open. 2024;6(8):951-958. doi:10.1016/j.cjco.2024.04.006

3. Anderson LN, Islam R, Sweetman A. Trends in obesity defined by body mass index among adults before and during the COVID-19 pandemic: a repeated cross-sectional study of the 2009-2023 Canadian Community Heath Surveys. CMAJ. 2025;197(25):E708-E718. doi:10.1503/cmaj.241421

4. Mo F, Pogany LM, Li FC, Morrison H. Prevalence of diabetes and cardiovascular comorbidity in the Canadian Community Health Survey 2002-2003. ScientificWorldJournal. 2006;6:96-105. doi:10.1100/tsw.2006.13

5. Geda NR, Janzen B, Pahwa P. Chronic disease multimorbidity among the Canadian population: prevalence and associated lifestyle factors. Arch Public Health. 2021;79(1):60. doi:10.1186/s13690-021-00583-7

6. Xiao X, Beach J, Senthilselvan A. Prevalence and determinants of multimorbidity in the Canadian population. PLoS One. 2024;19(1):e0297221. doi:10.1371/journal.pone.0297221

7. Azizi Z, Linder S, Ruiz MDCM, Alipour P, Raparelli V, Norris CM, et al. Assessing the relationship between sex, gender, and hypertension: a federated analysis of European and Canadian public health surveys. Medicine (Baltimore). 2025;104(6):e41021. doi:10.1097/MD.0000000000041021

8. Gisinger T, Azizi Z, Alipour P, Harreiter J, Raparelli V, Kublickiene K, et al. Sex and gender aspects in diabetes mellitus: focus on access to health care and cardiovascular outcomes. Front Public Health. 2023;11:1090541. doi:10.3389/fpubh.2023.1090541

9. Chukwudulue CR, Olawoye A, Igetei JE, Adeyanju FM, Oladimeji A, Arisoyin AE. Association between lifestyle factors and self-reported cardiovascular disease among Canadian adults aged 35 years and older: a study based on data from the Canadian Community Health Survey. Cureus. 2026;18(2):e104372. doi:10.7759/cureus.104372

10. Rothman KJ. The estimation of synergy or antagonism. Am J Epidemiol. 1976;103(5):506-511. doi:10.1093/oxfordjournals.aje.a112252

11. Knol MJ, VanderWeele TJ. Recommendations for presenting analyses of effect modification and interaction. Int J Epidemiol. 2012;41(2):514-520. doi:10.1093/ije/dyr218

12. De Angelis G, Chen Y. Obesity among women may increase the risk of arthritis: observations from the Canadian Community Health Survey, 2007-2008. Rheumatol Int. 2013;33(9):2249-2253. doi:10.1007/s00296-013-2712-5

13. St Sauver JL, Grossardt BR, Chamberlain AM, Kapoor E, Rocca WA. Synergistic interactions of obesity with sex, education, and smoking and accumulation of multi-morbidity (MM) across the lifespan. J Multimorb Comorb. 2023;13:26335565231160139. doi:10.1177/26335565231160139

14. Connor Gorber S, Shields M, Tremblay MS, McDowell I. The feasibility of establishing correction factors to adjust self-reported estimates of obesity. Health Rep. 2008;19(3):71-82.

15. Shields M, Connor Gorber S, Janssen I, Tremblay MS. Bias in self-reported estimates of obesity in Canadian health surveys: an update on correction equations for adults. Health Rep. 2011;22(3):35-45.

16. van Buuren S, Groothuis-Oudshoorn K. mice: multivariate imputation by chained equations in R. J Stat Softw. 2011;45(3). doi:10.18637/jss.v045.i03

17. Lumley T. Analysis of complex survey samples. J Stat Softw. 2004;9(8). doi:10.18637/jss.v009.i08

18. Li KH, Raghunathan TE, Rubin DB. Large-sample significance levels from multiply imputed data using moment-based statistics and an F reference distribution. J Am Stat Assoc. 1991;86(416):1065-1073. doi:10.1080/01621459.1991.10475152

19. Knol MJ, VanderWeele TJ, Groenwold RH, Klungel OH, Rovers MM, Grobbee DE. Estimating measures of interaction on an additive scale for preventive exposures. Eur J Epidemiol. 2011;26(6):433-438. doi:10.1007/s10654-011-9554-9

20. VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. Ann Intern Med. 2017;167(4):268-274. doi:10.7326/M16-2607

21. von Elm E, Altman DG, Egger M, Pocock SJ, Gotzsche PC, Vandenbroucke JP, et al. The Strengthening the Reporting of Observational Studies in Epidemiology (STROBE) statement: guidelines for reporting observational studies. Lancet. 2007;370(9596):1453-1457. doi:10.1016/S0140-6736(07)61602-X

22. Statistics Canada. Canadian Community Health Survey, Annual Component, 2017-2018 Public Use Microdata File (catalogue 82M0013X2024001). Ottawa: Statistics Canada; 2025.

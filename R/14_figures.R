# Figures for the report.
#   Figure 1: sex-specific odds ratios for overweight and obesity, all outcomes
#   Figure 2: standardized prevalence difference in multimorbidity associated
#             with overweight and obesity, by sex and age group
# Outputs: output/figures/figure1_forest.png, output/figures/figure2_additive.png

source("R/00_setup.R")
suppressPackageStartupMessages(library(ggplot2))

col_sex   <- c(Women = "#2a78d6", Men = "#eb6834")
shape_sex <- c(Women = 16, Men = 17)
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e6e5e1"

theme_fig <- function() {
  theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
          axis.text = element_text(colour = ink2), axis.title = element_text(colour = ink),
          strip.text = element_text(colour = ink, face = "bold", hjust = 0),
          legend.position = "top", legend.justification = "left",
          legend.title = element_blank(), legend.text = element_text(colour = ink),
          plot.caption = element_text(colour = ink2, hjust = 0, size = 8),
          plot.background = element_rect(fill = "white", colour = NA))
}

# Figure 1 ---------------------------------------------------------------------

outcomes <- c(cmm = "Multimorbidity", htn = "Hypertension", dm = "Diabetes", hd = "Heart disease")
ors <- fread(file.path(tab_dir, "models_interaction.csv"))[contrast %in% c("men", "women") & level %in% c("Overweight", "Obese")]
ratio <- fread(file.path(tab_dir, "models_interaction.csv"))[contrast == "ratio" & level %in% c("Overweight", "Obese")]
ors[, `:=`(sex = factor(fifelse(contrast == "women", "Women", "Men"), c("Women", "Men")),
           outcome_lab = factor(outcomes[outcome], outcomes),
           bmi = factor(level, c("Obese", "Overweight")))]
ratio[, `:=`(outcome_lab = factor(outcomes[outcome], outcomes), bmi = factor(level, c("Obese", "Overweight")),
             label = sprintf("Ratio %.2f (%.2f-%.2f)", or, lcl, ucl))]

p1 <- ggplot(ors, aes(x = or, y = bmi, colour = sex, shape = sex)) +
  geom_vline(xintercept = 1, colour = ink2, linewidth = 0.4) +
  geom_errorbar(aes(xmin = lcl, xmax = ucl), width = 0, linewidth = 0.6, orientation = "y",
                position = position_dodge(width = 0.55)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.55)) +
  geom_text(data = ratio, aes(x = 9.5, y = bmi, label = label), inherit.aes = FALSE,
            hjust = 1, size = 2.7, colour = ink2) +
  facet_wrap(~ outcome_lab, ncol = 1) +
  scale_x_log10(breaks = c(1, 1.5, 2, 3, 4, 6), limits = c(0.85, 9.5)) +
  scale_colour_manual(values = col_sex) + scale_shape_manual(values = shape_sex) +
  labs(x = "Odds ratio vs normal weight in the same sex (log scale)", y = NULL,
       caption = paste("Survey-weighted logistic regression with a sex x BMI interaction, adjusted for age group,",
                       "income,\neducation, smoking, alcohol and physical activity; 95% CIs from 1,000 bootstrap replicates",
                       "pooled over\n20 imputations. Ratio: odds ratio in women divided by odds ratio in men.")) +
  theme_fig() + theme(panel.grid.major.y = element_blank())
ggsave(file.path(fig_dir, "figure1_forest.png"), p1, width = 7, height = 7.2, dpi = 300, bg = "white")

# Figure 2 ---------------------------------------------------------------------

strata <- fread(file.path(tab_dir, "age_strata_additive.csv"))[outcome == "cmm" & measure %in% c("rd_female", "rd_male", "rd_diff")]
whole <- fread(file.path(tab_dir, "sensitivity.csv"))[variant == "excl_underweight" & measure %in% c("rd_female", "rd_male", "rd_diff")]
whole[, `:=`(stratum = "All ages", sparse = FALSE)]
fig2 <- rbind(strata[, .(stratum, bmi, measure, value, lcl, ucl, sparse)],
              whole[, .(stratum, bmi, measure, value, lcl, ucl, sparse)])
fig2[, `:=`(value = 100 * value, lcl = 100 * lcl, ucl = 100 * ucl,
            stratum = factor(stratum, c("18-44", "45-64", "65+", "All ages")),
            bmi = factor(bmi, c("Overweight", "Obese"), c("Overweight vs normal weight", "Obesity vs normal weight")))]
pts <- fig2[measure != "rd_diff"]
pts[, sex := factor(fifelse(measure == "rd_female", "Women", "Men"), c("Women", "Men"))]
diffs <- fig2[measure == "rd_diff"]
diffs[, label := sprintf("%+.1f\n(%+.1f to %+.1f)%s", value, lcl, ucl, fifelse(sparse, "*", ""))]
y_lab <- min(pts$lcl) - 3.2

p2 <- ggplot(pts, aes(x = stratum, y = value, colour = sex, shape = sex)) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = 0.4) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0, linewidth = 0.6, position = position_dodge(width = 0.5)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.5)) +
  geom_text(data = diffs, aes(x = stratum, y = y_lab, label = label), inherit.aes = FALSE,
            size = 2.5, colour = ink2, lineheight = 0.9) +
  facet_wrap(~ bmi, nrow = 1) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_y_continuous(breaks = seq(0, 20, 5)) +
  scale_colour_manual(values = col_sex) + scale_shape_manual(values = shape_sex) +
  labs(x = "Age group", y = "Difference in multimorbidity prevalence (percentage points)",
       caption = paste("Points: model-standardized prevalence difference against normal weight within each sex (95% CI).",
                       "\nNumbers below each group: women minus men (95% CI); negative values indicate a larger absolute excess",
                       "in men\n(additive interaction). Underweight excluded. * Fewer than 30 unweighted cases in a reference",
                       "cell; not interpreted.")) +
  theme_fig() + theme(panel.grid.major.x = element_blank())
ggsave(file.path(fig_dir, "figure2_additive.png"), p2, width = 7.5, height = 4.8, dpi = 300, bg = "white")
cat("Figures written\n")

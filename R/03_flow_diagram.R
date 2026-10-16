# Participant flow diagram from output/tables/flow_counts.csv.
# Output: output/figures/flow_diagram.png

source("R/00_setup.R")
suppressPackageStartupMessages(library(ggplot2))

flow <- fread(file.path(tab_dir, "flow_counts.csv"))
n <- setNames(flow$n, trimws(flow$step))
fmt <- function(x) format(x, big.mark = ",")

boxes <- data.table(
  x = c(0, 0, 0, 1.15, 1.15),
  y = c(4, 2.5, 0.5, 3.25, 1.5),
  w = c(1.5, 1.5, 1.5, 1.3, 1.3),
  h = c(0.7, 0.7, 0.7, 0.75, 1.1),
  label = c(
    sprintf("CCHS 2017-2018 PUMF respondents\nn = %s", fmt(n["All PUMF respondents"])),
    sprintf("Eligible cohort: adults 18+\nn = %s", fmt(n["Eligible cohort (adults 18+)"])),
    sprintf("Complete-case sample\nn = %s", fmt(n["Complete-case sample"])),
    sprintf("Excluded\nAged 12-17: %s\nCurrently pregnant: %s",
            fmt(n["Excluded: aged 12-17"]), fmt(n["Excluded: currently pregnant"])),
    sprintf("Incomplete records (retained\nfor multiple imputation)\nMissing outcome: %s\nMissing BMI, outcome known: %s\nMissing covariate only: %s",
            fmt(n["of which missing outcome"]), fmt(n["of which missing BMI category (outcome known)"]),
            fmt(n["of which missing covariate only"]))
  )
)

arrows <- data.table(
  x = c(0, 0, 0, 0), xend = c(0, 0, 0.5, 0.5),
  y = c(3.65, 2.15, 3.25, 1.5), yend = c(2.85, 0.85, 3.25, 1.5)
)

p <- ggplot() +
  geom_rect(data = boxes, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2),
            fill = "white", colour = "grey20", linewidth = 0.4) +
  geom_text(data = boxes, aes(x, y, label = label), size = 3.2, lineheight = 1.05) +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(2, "mm"), type = "closed"), linewidth = 0.4) +
  coord_cartesian(xlim = c(-0.8, 1.85), ylim = c(-0.1, 4.5)) +
  theme_void()

ggsave(file.path(fig_dir, "flow_diagram.png"), p, width = 7, height = 6, dpi = 300, bg = "white")

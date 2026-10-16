# Paths and constants shared by all scripts.
# Run scripts from the project root, e.g. `Rscript R/01_import.R`.

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

raw_dir     <- here("data", "raw", "2017-2018", "3226_CCHS-Ann_ESCC_2017-2018")
pumf_csv    <- file.path(raw_dir, "CCHS-82M0021X-E-2017-2018-Annual_F1.csv")
bsw_txt     <- file.path(raw_dir, "bsw", "Bootstrap weights_Poids bootstrap", "data", "bsw.txt")
derived_dir <- here("data", "derived")
out_dir     <- here("output")
fig_dir     <- here("output", "figures")
tab_dir     <- here("output", "tables")

for (d in c(derived_dir, out_dir, fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

n_boot <- 1000L

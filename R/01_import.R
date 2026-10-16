# Import the CCHS 2017-2018 PUMF and its bootstrap weights.
# Variable names and codes are from the 2017-2018 PUMF Data Dictionary.
# Outputs (not tracked): data/derived/pumf_selected.rds, data/derived/bsw.rds

source("R/00_setup.R")
suppressPackageStartupMessages(library(readr))

keep <- c(
  "ADM_RNO", "GEO_PRV", "ADM_PRX", "DHH_SEX", "DHHGAGE",
  "HWTDGBMI", "HWTDGCOR", "HWTDGISW", "HWTDGBCC",
  "CCC_065", "CCC_070", "CCC_085", "CCC_090", "CCC_095", "CCC_105", "CCC_110",
  "EHG2DVH3", "INCDGHH", "SMKDVSTY", "ALCDVTTM", "PAADVWHO", "PAADVVOL",
  "WTS_M"
)

pumf <- fread(pumf_csv, select = keep)
stopifnot(nrow(pumf) == 113290L, !anyDuplicated(pumf$ADM_RNO))
saveRDS(pumf, file.path(derived_dir, "pumf_selected.rds"))

# Bootstrap file layout (BSW_i.sas in the release package):
# ADM_RNO cols 1-20, FWGT cols 21-27 (7.2), BSW1-BSW1000 7.2 each from col 28.
# No decimal points appear in the file, so the two decimals are implied.
bsw <- read_fwf(
  bsw_txt,
  col_positions = fwf_widths(c(20, 7, rep(7, n_boot)), c("ADM_RNO", "FWGT", paste0("BSW", seq_len(n_boot)))),
  col_types = cols(.default = col_integer()),
  progress = FALSE
)
setDT(bsw)
stopifnot(nrow(bsw) == nrow(pumf), !anyDuplicated(bsw$ADM_RNO))

wcols <- c("FWGT", paste0("BSW", seq_len(n_boot)))
bsw[, (wcols) := lapply(.SD, function(x) x / 100), .SDcols = wcols]

chk <- merge(pumf[, .(ADM_RNO, WTS_M)], bsw[, .(ADM_RNO, FWGT)], by = "ADM_RNO")
stopifnot(nrow(chk) == nrow(pumf), max(abs(chk$WTS_M - chk$FWGT)) < 0.01)
stopifnot(min(as.matrix(bsw[, ..wcols])) > 0)

saveRDS(bsw, file.path(derived_dir, "bsw.rds"))

cat("PUMF records:", nrow(pumf), "\n")
cat("Bootstrap records:", nrow(bsw), " replicates:", n_boot, "\n")
cat("Max |WTS_M - FWGT|:", max(abs(chk$WTS_M - chk$FWGT)), "\n")

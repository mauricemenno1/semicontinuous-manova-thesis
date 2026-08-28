###############################################################################
## Final RNA-calibrated simulation settings
###############################################################################
source(file.path("R", "simulation_common.R"))
source(file.path("R", "fit_scmanova.R"))

B <- 1000L
NREP <- 1000L
K <- 2L

## The thesis specifies the same search interval for both occurrence models.
LAMBDA_MV <- c(0, 100)
LAMBDA0_MV <- c(0, 100)
LAMBDA_BERNOULLI <- c(0, 100)
LAMBDA0_BERNOULLI <- c(0, 100)

STUDY_DIR <- file.path("simulations", "rna_calibrated")
RAW_DIR <- file.path(STUDY_DIR, "results", "raw")
SUMMARY_FILE <- file.path(STUDY_DIR, "results", "summary.csv")
TABLE_DIR <- file.path(STUDY_DIR, "tables")
DATA_FILE <- file.path("data", "rna_simulation_parameters.RData")

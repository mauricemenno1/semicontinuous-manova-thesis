###############################################################################
## Final thesis simulation: exchangeable multivariate-Bernoulli injection
###############################################################################
source(file.path("R", "simulation_common.R"))
source(file.path("R", "fit_scmanova.R"))
source(file.path("R", "fit_chen.R"))

K <- 2L
B <- 100L
NREP <- 250L
P_VALUES <- c(50L, 100L, 150L)
N_PER_GROUP_VALUES <- c(5L, 10L)
PMISS_VALUES <- c(0.2, 0.5, 0.8)
RHO_VALUES <- c(0, 0.4)
OCCURRENCE_ICC <- 0.05

STUDY_DIR <- file.path("simulations", "exchangeable_occurrence")
RAW_DIR <- file.path(STUDY_DIR, "results", "raw")
SUMMARY_FILE <- file.path(STUDY_DIR, "results", "summary.csv")
TABLE_DIR <- file.path(STUDY_DIR, "tables")

GRID <- build_scenario_grid(
  p_values = P_VALUES,
  n_per_group_values = N_PER_GROUP_VALUES,
  pmiss_values = PMISS_VALUES,
  rho_values = RHO_VALUES,
  K = K
)

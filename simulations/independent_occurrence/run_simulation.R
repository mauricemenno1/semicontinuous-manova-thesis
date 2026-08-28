###############################################################################
## Run all independent-Bernoulli-injection scenarios serially on one core
###############################################################################
source(file.path("simulations", "independent_occurrence", "config.R"))
assert_repository_root()
assert_modified_package(require_bernoulli = TRUE)

independent_data_generator <- function(scenario, replication_index) {
  set.seed(1235L + as.integer(replication_index))
  pars <- scenario_parameters(scenario, K = K)

  x <- semicontMANOVA::scMANOVAsimulation(
    n = pars$n,
    p = pars$p,
    pmiss = pars$pmiss,
    rho = as.numeric(scenario$rho),
    mu = pars$mu
  )

  ## Restore this exact post-simulation RNG state before each model fit. This
  ## mirrors the original replication scripts while ensuring every method gets
  ## the same generated data and the same starting resampling stream.
  list(
    x = x,
    n = pars$n,
    fit_rng = list(state = .Random.seed)
  )
}

MODEL_SPECS <- list(
  CHEN = list(
    fit = function(x, n, B, fit_rng) {
      fit_chen_for_study(x, n, B = B, fit_rng = fit_rng)
    }
  ),
  scMAN_MV = list(
    fit = function(x, n, B, fit_rng) {
      fit_scmanova_for_study(
        x, n, B = B, missing_model = NULL,
        lambda = c(0, 100), lambda0 = c(0, 100), fit_rng = fit_rng
      )
    }
  ),
  scMAN_Bernoulli = list(
    fit = function(x, n, B, fit_rng) {
      fit_scmanova_for_study(
        x, n, B = B, missing_model = "bernoulli",
        lambda = c(0, 100), lambda0 = c(0, 100), fit_rng = fit_rng
      )
    }
  )
)

cat("Independent Bernoulli injection study\n")
cat("Serial execution: 1 core\n")
cat("B =", B, "| nrep =", NREP, "| p =", paste(P_VALUES, collapse = ", "), "\n")

run_serial_study(
  grid = GRID,
  nrep = NREP,
  B = B,
  raw_dir = RAW_DIR,
  model_specs = MODEL_SPECS,
  data_generator = independent_data_generator,
  progress_every = 25L
)

###############################################################################
## Run the RNA-calibrated H0/H1 simulations serially on one core
###############################################################################
source(file.path("simulations", "rna_calibrated", "config.R"))
assert_repository_root()
assert_modified_package(require_bernoulli = TRUE)

if (!file.exists(DATA_FILE)) stop("Missing data file: ", DATA_FILE)

e <- new.env(parent = emptyenv())
load(DATA_FILE, envir = e)

## Apply the same preprocessing as the supplied replication code.
for (k in seq_along(e$n)) {
  e$mu1[k, which(is.na(e$mu1[k, ]))] <- e$mu0[which(is.na(e$mu1[k, ]))]
}

e$sigma0[which(is.na(diag(e$sigma0))), which(is.na(diag(e$sigma0)))] <- 1
e$sigma0[which(is.na(e$sigma0))] <- 0
diag(e$sigma0) <- diag(e$sigma0) + 0.1

e$sigma1[which(is.na(diag(e$sigma1))), which(is.na(diag(e$sigma1)))] <- 1
e$sigma1[which(is.na(e$sigma1))] <- 0
diag(e$sigma1) <- diag(e$sigma1) + 0.1

SCENARIOS <- list(
  H0 = list(mu = e$mu0, sigma = e$sigma0, pmiss = rep(e$pmiss0, K)),
  H1 = list(mu = e$mu1, sigma = e$sigma1, pmiss = e$pmiss1)
)

MODEL_SPECS <- list(
  scMAN_MV = list(
    fit = function(x, n, B, fit_rng) {
      fit_scmanova_for_study(
        x, n, B = B, missing_model = NULL,
        lambda = LAMBDA_MV, lambda0 = LAMBDA0_MV, fit_rng = fit_rng
      )
    }
  ),
  scMAN_Bernoulli = list(
    fit = function(x, n, B, fit_rng) {
      fit_scmanova_for_study(
        x, n, B = B, missing_model = "bernoulli",
        lambda = LAMBDA_BERNOULLI,
        lambda0 = LAMBDA0_BERNOULLI,
        fit_rng = fit_rng
      )
    }
  )
)

dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
cat("RNA-calibrated simulation\n")
cat("Serial execution: 1 core\n")
cat("B =", B, "| nrep =", NREP, "| n =", paste(e$n, collapse = "+"), "| p =", e$p, "\n")

for (scenario_name in names(SCENARIOS)) {
  raw_path <- file.path(RAW_DIR, paste0(scenario_name, ".csv"))
  results <- read_raw_results(raw_path)
  pars <- SCENARIOS[[scenario_name]]

  cat("\nScenario ", scenario_name, " | existing rows: ", nrow(results), "\n", sep = "")

  for (irep in seq_len(NREP)) {
    done_models <- results$model[results$replication == irep]
    missing_models <- setdiff(names(MODEL_SPECS), done_models)
    if (length(missing_models) == 0L) next

    set.seed(1235L + irep)
    x <- semicontMANOVA::scMANOVAsimulation(
      n = e$n,
      p = e$p,
      pmiss = pars$pmiss,
      sigma = pars$sigma,
      mu = pars$mu
    )
    fit_rng <- list(state = .Random.seed)

    for (model_name in missing_models) {
      fit <- MODEL_SPECS[[model_name]]$fit(
        x = x,
        n = e$n,
        B = B,
        fit_rng = fit_rng
      )
      row <- cbind(
        data.frame(replication = irep, model = model_name, stringsAsFactors = FALSE),
        fit,
        data.frame(retained_p = ncol(x))
      )
      results <- rbind(results, row)
    }

    results <- results[order(results$replication, results$model), , drop = FALSE]
    rownames(results) <- NULL
    write_csv_atomic(results, raw_path)

    if (irep %% 25L == 0L || irep == NREP) {
      cat("  replication ", irep, "/", NREP, " complete\n", sep = "")
      flush.console()
    }
  }
}

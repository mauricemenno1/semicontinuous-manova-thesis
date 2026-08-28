###############################################################################
## Run all exchangeable-Bernoulli-injection scenarios serially on one core
###############################################################################
source(file.path("simulations", "exchangeable_occurrence", "config.R"))
assert_repository_root()
assert_modified_package(require_bernoulli = TRUE)

if (!requireNamespace("mvtnorm", quietly = TRUE)) {
  stop("Package 'mvtnorm' is required.")
}

positive_count_distribution <- function(p, pmiss, icc) {
  if (pmiss <= 0) {
    out <- numeric(p + 1L); out[p + 1L] <- 1; return(out)
  }
  if (pmiss >= 1) {
    out <- numeric(p + 1L); out[1L] <- 1; return(out)
  }

  concentration <- (1 / icc) - 1
  shape_positive <- (1 - pmiss) * concentration
  shape_zero <- pmiss * concentration
  s <- 0:p

  log_probability <-
    lchoose(p, s) +
    lbeta(s + shape_positive, p - s + shape_zero) -
    lbeta(shape_positive, shape_zero)

  log_probability <- log_probability - max(log_probability)
  probability <- exp(log_probability)
  probability / sum(probability)
}

simulate_exchangeable_presence <- function(n_obs, p, pmiss, icc) {
  count_probability <- positive_count_distribution(p, pmiss, icc)
  positive_counts <- sample.int(
    p + 1L,
    size = n_obs,
    replace = TRUE,
    prob = count_probability
  ) - 1L

  y <- matrix(0L, nrow = n_obs, ncol = p)
  for (i in seq_len(n_obs)) {
    s <- positive_counts[i]
    if (s > 0L) {
      y[i, sample.int(p, size = s, replace = FALSE)] <- 1L
    }
  }
  y
}

simulate_exchangeable_dataset <- function(n, p, pmiss, rho, mu, icc) {
  K_local <- length(n)
  if (is.null(mu)) mu <- matrix(0, nrow = K_local, ncol = p)
  if (length(pmiss) != K_local) pmiss <- rep(pmiss, length.out = K_local)

  sigma <- matrix(rho, nrow = p, ncol = p)
  diag(sigma) <- 1

  original <- matrix(NA_real_, nrow = sum(n), ncol = p)
  y <- matrix(1L, nrow = sum(n), ncol = p)
  start <- 1L

  for (k in seq_len(K_local)) {
    end <- start + n[k] - 1L
    original[start:end, ] <- mvtnorm::rmvnorm(
      n[k],
      mean = mu[k, ],
      sigma = sigma
    )
    y[start:end, ] <- simulate_exchangeable_presence(
      n_obs = n[k], p = p, pmiss = pmiss[k], icc = icc
    )
    start <- end + 1L
  }

  x <- original
  x[y == 0L] <- 0

  ## Match scMANOVAsimulation(): discard components that are zero everywhere.
  retained <- which(colSums(y) != 0L)
  x[, retained, drop = FALSE]
}

exchangeable_data_generator <- function(scenario, replication_index) {
  set.seed(1235L + as.integer(replication_index))
  pars <- scenario_parameters(scenario, K = K)
  x <- simulate_exchangeable_dataset(
    n = pars$n,
    p = pars$p,
    pmiss = pars$pmiss,
    rho = as.numeric(scenario$rho),
    mu = pars$mu,
    icc = OCCURRENCE_ICC
  )

  ## This is the seed convention used in the final multivariate-injection run.
  fit_seed <- as.integer(
    100000L + 1000L * as.integer(scenario$scenario_id) +
      as.integer(replication_index)
  )

  list(x = x, n = pars$n, fit_rng = list(seed = fit_seed))
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

cat("Exchangeable multivariate-Bernoulli injection study\n")
cat("Serial execution: 1 core\n")
cat("Occurrence ICC =", OCCURRENCE_ICC, "\n")
cat("B =", B, "| nrep =", NREP, "| p =", paste(P_VALUES, collapse = ", "), "\n")

run_serial_study(
  grid = GRID,
  nrep = NREP,
  B = B,
  raw_dir = RAW_DIR,
  model_specs = MODEL_SPECS,
  data_generator = exchangeable_data_generator,
  progress_every = 25L
)

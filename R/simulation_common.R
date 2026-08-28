###############################################################################
## Shared utilities for the thesis simulation studies
###############################################################################

assert_repository_root <- function() {
  required <- c(
    file.path("R", "fit_scmanova.R"),
    file.path("R", "fit_chen.R"),
    file.path("package", "semicontMANOVA_0.2_bernoulli.tar.gz")
  )
  missing <- required[!file.exists(required)]
  if (length(missing) > 0L) {
    stop(
      "Run this script from the repository root. Missing: ",
      paste(missing, collapse = ", ")
    )
  }
}

assert_modified_package <- function(require_bernoulli = TRUE) {
  if (!requireNamespace("semicontMANOVA", quietly = TRUE)) {
    stop(
      "Package 'semicontMANOVA' is not installed. Run ",
      "source('package/install_package.R') first."
    )
  }
  if (isTRUE(require_bernoulli)) {
    has_arg <- "missing.model" %in% names(formals(semicontMANOVA::scMANOVA))
    if (!has_arg) {
      stop(
        "The installed semicontMANOVA package is not the edited thesis version: ",
        "scMANOVA() has no 'missing.model' argument."
      )
    }
  }
}

scenario_types <- list(
  H0 = c(c1 = 0, c2 = 0),
  H1_mean_1 = c(c1 = 1, c2 = 0),
  H1_mean_5 = c(c1 = 5, c2 = 0),
  H1_missing_0.15 = c(c1 = 0, c2 = 0.15),
  H1_missing_0.30 = c(c1 = 0, c2 = 0.30)
)

build_scenario_grid <- function(
    p_values = c(50L, 100L, 150L),
    n_per_group_values = c(5L, 10L),
    pmiss_values = c(0.2, 0.5, 0.8),
    rho_values = c(0, 0.4),
    K = 2L
) {
  rows <- list()
  id <- 0L

  ## p is the outer loop deliberately: a serial run therefore completes all
  ## p = 50 scenarios before p = 100, and p = 100 before p = 150.
  for (p in p_values) {
    for (type_name in names(scenario_types)) {
      c1 <- unname(scenario_types[[type_name]]["c1"])
      c2 <- unname(scenario_types[[type_name]]["c2"])

      for (n_per_group in n_per_group_values) {
        for (pmiss in pmiss_values) {
          if ((pmiss + c2) >= 1) next
          for (rho in rho_values) {
            id <- id + 1L
            rows[[id]] <- data.frame(
              scenario_id = id,
              scenario = type_name,
              c1 = c1,
              c2 = c2,
              n_per_group = as.integer(n_per_group),
              n_total = as.integer(K * n_per_group),
              p = as.integer(p),
              pmiss = as.numeric(pmiss),
              rho = as.numeric(rho),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

scenario_key <- function(scenario) {
  paste0(
    scenario$scenario,
    "__n", scenario$n_total,
    "__p", scenario$p,
    "__pmiss", format(scenario$pmiss, trim = TRUE, scientific = FALSE),
    "__rho", format(scenario$rho, trim = TRUE, scientific = FALSE)
  )
}

scenario_parameters <- function(scenario, K = 2L) {
  n <- rep(as.integer(scenario$n_per_group), K)
  p <- as.integer(scenario$p)
  c1 <- as.numeric(scenario$c1)
  c2 <- as.numeric(scenario$c2)
  base_pmiss <- as.numeric(scenario$pmiss)

  if (c1 == 0 && c2 == 0) {
    mu <- NULL
    pmiss <- rep(base_pmiss, K)
  } else if (c1 > 0 && c2 == 0) {
    group_shift <- c1 * seq(0, K - 1) / (K - 1)
    mu <- matrix(group_shift, nrow = K, ncol = p)
    pmiss <- rep(base_pmiss, K)
  } else if (c1 == 0 && c2 > 0) {
    mu <- NULL
    pmiss <- base_pmiss + c2 * seq(0, K - 1) / (K - 1)
  } else {
    stop("Unsupported c1/c2 combination.")
  }

  list(n = n, p = p, mu = mu, pmiss = pmiss)
}

empty_result_row <- function() {
  data.frame(
    replication = integer(),
    model = character(),
    elapsed_seconds = numeric(),
    logLik = numeric(),
    logLik0 = numeric(),
    lambda = numeric(),
    lambda0 = numeric(),
    df = numeric(),
    df0 = numeric(),
    aic = numeric(),
    aic0 = numeric(),
    statistic = numeric(),
    p_value = numeric(),
    dimSigma = numeric(),
    dimSigma0 = numeric(),
    retained_p = integer(),
    stringsAsFactors = FALSE
  )
}

read_raw_results <- function(path) {
  if (!file.exists(path) || isTRUE(file.info(path)$size == 0)) {
    return(empty_result_row())
  }
  x <- read.csv(path, stringsAsFactors = FALSE)
  expected <- names(empty_result_row())
  missing <- setdiff(expected, names(x))
  if (length(missing) > 0L) {
    stop("Raw result file has unexpected columns: ", path)
  }
  x[, expected, drop = FALSE]
}

write_csv_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  write.csv(x, tmp, row.names = FALSE, na = "NA")
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) {
    ok <- file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
    if (!isTRUE(ok)) stop("Could not write: ", path)
  }
  invisible(path)
}

set_fit_rng <- function(fit_rng) {
  if (is.null(fit_rng)) return(invisible(NULL))
  if (!is.null(fit_rng$state)) {
    assign(".Random.seed", fit_rng$state, envir = .GlobalEnv)
  } else if (!is.null(fit_rng$seed)) {
    set.seed(as.integer(fit_rng$seed))
  } else {
    stop("fit_rng must contain either $state or $seed.")
  }
  invisible(NULL)
}

fit_scmanova_for_study <- function(
    x,
    n,
    B,
    missing_model,
    lambda = c(0, 100),
    lambda0 = c(0, 100),
    fit_rng = NULL
) {
  set_fit_rng(fit_rng)
  timing <- system.time({
    value <- estimation(
      x = x,
      n = n,
      lambda = lambda,
      lambda0 = lambda0,
      B = B,
      penalty = function(n, p) log(n) + 0.5 * log(p),
      ident = TRUE,
      missing.model = missing_model
    )
  })

  value <- as.numeric(value)
  if (length(value) < 12L) value <- c(value, rep(NA_real_, 12L - length(value)))
  value <- value[seq_len(12L)]

  data.frame(
    elapsed_seconds = unname(timing[["elapsed"]]),
    logLik = value[1], logLik0 = value[2],
    lambda = value[3], lambda0 = value[4],
    df = value[5], df0 = value[6],
    aic = value[7], aic0 = value[8],
    statistic = value[9], p_value = value[10],
    dimSigma = value[11], dimSigma0 = value[12]
  )
}

fit_chen_for_study <- function(x, n, B, fit_rng = NULL) {
  set_fit_rng(fit_rng)
  timing <- system.time({
    value <- try(
      estimationCHEN(
        x = x,
        n = n,
        Lambda = seq(0.01, 10, 0.05),
        B = B,
        delta = 0.5,
        alpha = 0.05
      ),
      silent = TRUE
    )
  })
  if (inherits(value, "try-error") || length(value) != 1L || !is.finite(value)) {
    value <- NA_real_
  }

  data.frame(
    elapsed_seconds = unname(timing[["elapsed"]]),
    logLik = NA_real_, logLik0 = NA_real_,
    lambda = NA_real_, lambda0 = NA_real_,
    df = NA_real_, df0 = NA_real_,
    aic = NA_real_, aic0 = NA_real_,
    statistic = NA_real_, p_value = as.numeric(value),
    dimSigma = NA_real_, dimSigma0 = NA_real_
  )
}

run_serial_study <- function(
    grid,
    nrep,
    B,
    raw_dir,
    model_specs,
    data_generator,
    progress_every = 25L
) {
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  for (s in seq_len(nrow(grid))) {
    scenario <- grid[s, , drop = FALSE]
    key <- scenario_key(scenario)
    raw_path <- file.path(raw_dir, paste0(key, ".csv"))
    results <- read_raw_results(raw_path)

    cat(
      "\n[", s, "/", nrow(grid), "] ", key,
      " | existing rows: ", nrow(results), "\n",
      sep = ""
    )

    for (irep in seq_len(nrep)) {
      done_models <- results$model[results$replication == irep]
      missing_models <- setdiff(names(model_specs), done_models)
      if (length(missing_models) == 0L) next

      generated <- data_generator(scenario, irep)
      x <- generated$x
      n <- generated$n

      for (model_name in missing_models) {
        fit <- model_specs[[model_name]]$fit(
          x = x,
          n = n,
          B = B,
          fit_rng = generated$fit_rng
        )

        row <- cbind(
          data.frame(
            replication = irep,
            model = model_name,
            stringsAsFactors = FALSE
          ),
          fit,
          data.frame(retained_p = ncol(x))
        )
        results <- rbind(results, row)
      }

      results <- results[order(results$replication, results$model), , drop = FALSE]
      rownames(results) <- NULL
      write_csv_atomic(results, raw_path)

      if (irep %% progress_every == 0L || irep == nrep) {
        cat("  replication ", irep, "/", nrep, " complete\n", sep = "")
        flush.console()
      }
    }
  }

  invisible(TRUE)
}

summarise_study <- function(grid, raw_dir, output_file, alpha = 0.05) {
  rows <- list()

  for (s in seq_len(nrow(grid))) {
    scenario <- grid[s, , drop = FALSE]
    path <- file.path(raw_dir, paste0(scenario_key(scenario), ".csv"))
    x <- read_raw_results(path)
    if (nrow(x) == 0L) next

    for (model_name in unique(x$model)) {
      z <- x[x$model == model_name, , drop = FALSE]
      ok <- !is.na(z$p_value)
      rows[[length(rows) + 1L]] <- data.frame(
        scenario_id = scenario$scenario_id,
        scenario = scenario$scenario,
        c1 = scenario$c1,
        c2 = scenario$c2,
        n_per_group = scenario$n_per_group,
        n_total = scenario$n_total,
        p = scenario$p,
        pmiss = scenario$pmiss,
        rho = scenario$rho,
        model = model_name,
        completed = sum(ok),
        rejection_proportion = if (any(ok)) mean(z$p_value[ok] < alpha) else NA_real_,
        mean_lambda = if (all(is.na(z$lambda))) NA_real_ else mean(z$lambda, na.rm = TRUE),
        mean_lambda0 = if (all(is.na(z$lambda0))) NA_real_ else mean(z$lambda0, na.rm = TRUE),
        mean_dim = if (all(is.na(z$dimSigma))) NA_real_ else mean(z$dimSigma, na.rm = TRUE),
        mean_dim0 = if (all(is.na(z$dimSigma0))) NA_real_ else mean(z$dimSigma0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0L) {
    stop("No raw result files were found in: ", raw_dir)
  }

  out <- do.call(rbind, rows)
  out <- out[order(out$p, out$scenario_id, out$model), , drop = FALSE]
  rownames(out) <- NULL
  write_csv_atomic(out, output_file)
  out
}

format_table_value <- function(x) {
  if (length(x) == 0L || is.na(x)) return("x")
  sprintf("%.3f", x)
}

escape_latex <- function(x) {
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x
}

###############################################################################
## Helpers for rebuilding thesis tables from the archived historical outputs
###############################################################################

historical_scenario_map <- data.frame(
  scenario = c("h0", "mean1", "mean5", "prob015", "prob030"),
  c1 = c(0, 1, 5, 0, 0),
  c2 = c(0, 0, 0, 0.15, 0.30),
  stringsAsFactors = FALSE
)

historical_number_string <- function(x) {
  format(as.numeric(x), trim = TRUE, scientific = FALSE, digits = 15)
}

historical_raw_rejection <- function(path, model, alpha = 0.05) {
  x <- read.table(
    path,
    sep = ";",
    header = FALSE,
    fill = TRUE,
    blank.lines.skip = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(x) == 0L) stop("Historical raw result file is empty: ", path)
  if (ncol(x) == 1L) x <- t(x)

  p_col <- if (identical(model, "CHEN")) 2L else 16L
  if (ncol(x) < p_col) {
    stop("Historical raw result has too few columns: ", path)
  }

  p_values <- suppressWarnings(as.numeric(x[, p_col]))
  p_values <- p_values[!is.na(p_values)]
  if (length(p_values) == 0L) {
    stop("No completed permutation p-values found in: ", path)
  }

  list(
    rejection = round(mean(p_values < alpha), 3),
    n_completed = length(p_values)
  )
}

historical_find_raw_file <- function(method_dir, scenario, n_total, p, pmiss, rho) {
  raw_dir <- file.path(method_dir, scenario)
  if (!dir.exists(raw_dir)) return(NULL)

  filename <- paste0(
    "n", historical_number_string(n_total),
    "_p", historical_number_string(p),
    "_pi", historical_number_string(pmiss),
    "_rho", historical_number_string(rho),
    ".txt"
  )

  path <- file.path(raw_dir, filename)
  if (!file.exists(path)) return(NULL)
  path
}

collect_historical_simulation_summary <- function(historical_root, alpha = 0.05) {
  method_folders <- c(
    scMAN_MV = "scman_exc",
    scMAN_Bernoulli = "scman_ind",
    CHEN = "chen"
  )

  rows <- list()

  for (model in names(method_folders)) {
    method_dir <- file.path(historical_root, method_folders[[model]])
    if (!dir.exists(method_dir)) {
      stop("Historical method folder not found: ", method_dir)
    }

    for (i in seq_len(nrow(historical_scenario_map))) {
      scenario <- historical_scenario_map$scenario[i]
      summary_path <- file.path(method_dir, paste0(scenario, ".csv"))
      if (!file.exists(summary_path)) {
        stop("Historical summary file not found: ", summary_path)
      }

      x <- read.table(
        summary_path,
        sep = ";",
        header = TRUE,
        stringsAsFactors = FALSE,
        na.strings = c("NA", "")
      )

      ## Thesis grid: p = 50, 100, 150 and total n = 10, 20
      ## (n_k = 5, 10 in each of the two groups).
      x <- x[x$p %in% c(50, 100, 150) & x$n %in% c(10, 20), , drop = FALSE]

      for (j in seq_len(nrow(x))) {
        rejection <- suppressWarnings(as.numeric(x$rejPERM[j]))
        n_completed <- if ("n_simEnded" %in% names(x)) as.integer(x$n_simEnded[j]) else NA_integer_
        source_used <- "summary_csv"

        ## One archived independent-occurrence summary row is stale/NA even
        ## though its raw 250-replication file is complete. Keep the archived
        ## CSV unchanged and use the raw file only as a fallback.
        if (is.na(rejection)) {
          raw_path <- historical_find_raw_file(
            method_dir = method_dir,
            scenario = scenario,
            n_total = x$n[j],
            p = x$p[j],
            pmiss = x$pmiss[j],
            rho = x$rho[j]
          )
          if (is.null(raw_path)) {
            stop(
              "Historical summary is missing a rejection proportion and no raw ",
              "result file was found for ", model, ", ", scenario,
              ", n=", x$n[j], ", p=", x$p[j],
              ", pmiss=", x$pmiss[j], ", rho=", x$rho[j]
            )
          }
          raw_summary <- historical_raw_rejection(raw_path, model, alpha = alpha)
          rejection <- raw_summary$rejection
          n_completed <- raw_summary$n_completed
          source_used <- "raw_txt_fallback"
        }

        rows[[length(rows) + 1L]] <- data.frame(
          c1 = historical_scenario_map$c1[i],
          c2 = historical_scenario_map$c2[i],
          pmiss = as.numeric(x$pmiss[j]),
          rho = as.numeric(x$rho[j]),
          model = model,
          n_per_group = as.integer(x$n[j] / 2),
          p = as.integer(x$p[j]),
          rejection_proportion = rejection,
          historical_n_completed = n_completed,
          historical_source = source_used,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

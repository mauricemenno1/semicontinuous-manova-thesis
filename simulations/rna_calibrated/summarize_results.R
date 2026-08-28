###############################################################################
## Summarise the RNA-calibrated simulation
###############################################################################
source(file.path("simulations", "rna_calibrated", "config.R"))
assert_repository_root()

if (!file.exists(DATA_FILE)) stop("Missing data file: ", DATA_FILE)
e <- new.env(parent = emptyenv())
load(DATA_FILE, envir = e)

rows <- list()
for (scenario_name in c("H0", "H1")) {
  path <- file.path(RAW_DIR, paste0(scenario_name, ".csv"))
  x <- read_raw_results(path)
  if (nrow(x) == 0L) next

  for (model_name in unique(x$model)) {
    z <- x[x$model == model_name, , drop = FALSE]
    ok <- !is.na(z$p_value)
    rows[[length(rows) + 1L]] <- data.frame(
      model = model_name,
      scenario = scenario_name,
      n = sum(e$n),
      p = e$p,
      completed = sum(ok),
      rejection_proportion = if (any(ok)) mean(z$p_value[ok] < 0.05) else NA_real_,
      mean_dim = if (all(is.na(z$dimSigma))) NA_real_ else mean(z$dimSigma, na.rm = TRUE),
      mean_dim0 = if (all(is.na(z$dimSigma0))) NA_real_ else mean(z$dimSigma0, na.rm = TRUE),
      mean_lambda = if (all(is.na(z$lambda))) NA_real_ else mean(z$lambda, na.rm = TRUE),
      mean_lambda0 = if (all(is.na(z$lambda0))) NA_real_ else mean(z$lambda0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}

if (length(rows) == 0L) stop("No RNA raw results were found.")
summary_df <- do.call(rbind, rows)
summary_df <- summary_df[order(summary_df$model, summary_df$scenario), , drop = FALSE]
rownames(summary_df) <- NULL
write_csv_atomic(summary_df, SUMMARY_FILE)
cat("Wrote:", SUMMARY_FILE, "\n")
print(summary_df)

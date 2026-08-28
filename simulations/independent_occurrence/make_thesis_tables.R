###############################################################################
## Recreate the three submitted-thesis tables from archived historical results
###############################################################################
source(file.path("R", "historical_results.R"))
source(file.path("R", "make_simulation_tables.R"))

historical_root <- file.path("simulations", "independent_occurrence", "history")
out_dir <- file.path("simulations", "independent_occurrence", "tables", "thesis")

summary_df <- collect_historical_simulation_summary(historical_root)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  summary_df,
  file.path(out_dir, "source_summary.csv"),
  row.names = FALSE,
  na = "NA"
)
write_simulation_tables(summary_df, out_dir, study_key = "independent")

cat("Wrote independent-Bernoulli thesis tables to:", out_dir, "\n")

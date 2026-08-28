###############################################################################
## Summarise the independent-Bernoulli-injection simulation
###############################################################################
source(file.path("simulations", "independent_occurrence", "config.R"))
assert_repository_root()

summary_df <- summarise_study(
  grid = GRID,
  raw_dir = RAW_DIR,
  output_file = SUMMARY_FILE,
  alpha = 0.05
)
cat("Wrote:", SUMMARY_FILE, "\n")
print(summary_df)

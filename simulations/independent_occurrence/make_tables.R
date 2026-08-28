###############################################################################
## Make the three thesis tables for independent Bernoulli occurrence injection
###############################################################################
source(file.path("simulations", "independent_occurrence", "config.R"))
source(file.path("R", "make_simulation_tables.R"))
assert_repository_root()

if (!file.exists(SUMMARY_FILE)) stop("Run summarize_results.R first.")
summary_df <- read.csv(SUMMARY_FILE, stringsAsFactors = FALSE)
write_simulation_tables(summary_df, TABLE_DIR, study_key = "independent")
cat("Wrote thesis tables to:", TABLE_DIR, "\n")

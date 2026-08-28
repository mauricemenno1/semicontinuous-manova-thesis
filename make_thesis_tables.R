###############################################################################
## Recreate all seven submitted-thesis tables from archived result files
##
## Run from the repository root. No Monte Carlo simulations are rerun.
###############################################################################

source(file.path("simulations", "independent_occurrence", "make_thesis_tables.R"))
source(file.path("simulations", "exchangeable_occurrence", "make_thesis_tables.R"))
source(file.path("simulations", "rna_calibrated", "make_thesis_table.R"))

cat("\nAll seven thesis tables have been written from archived results.\n")

###############################################################################
## Run the complete Zeller colorectal-cancer application
## Run from the repository root.
###############################################################################

source(file.path("applications", "zeller_crc", "prepare_data.R"))
source(file.path("applications", "zeller_crc", "summarize_data.R"))
source(file.path("applications", "zeller_crc", "run_models.R"))
source(file.path("applications", "zeller_crc", "make_figures.R"))

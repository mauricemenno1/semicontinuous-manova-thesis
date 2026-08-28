###############################################################################
## Shared helpers for the Zeller colorectal-cancer application
###############################################################################

zeller_dir <- file.path("applications", "zeller_crc")
zeller_results_dir <- file.path(zeller_dir, "results")
zeller_figures_dir <- file.path(zeller_dir, "figures")

check_zeller_root <- function() {
  if (!file.exists(file.path(zeller_dir, "config.R"))) {
    stop("Run the Zeller scripts from the repository root.")
  }
  dir.create(zeller_results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(zeller_figures_dir, recursive = TRUE, showWarnings = FALSE)
}

write_csv <- function(x, filename) {
  utils::write.csv(
    x,
    file.path(zeller_results_dir, filename),
    row.names = FALSE
  )
}


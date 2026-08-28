###############################################################################
## Install and verify the edited semicontMANOVA package used in the thesis
###############################################################################

package_file <- file.path("package", "semicontMANOVA_0.2_bernoulli.tar.gz")
if (!file.exists(package_file)) {
  stop("Run this script from the repository root.")
}

needed <- c("matrixcalc", "mvtnorm")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  install.packages(missing)
}

install.packages(package_file, repos = NULL, type = "source")

has_missing_model <- "missing.model" %in%
  names(formals(semicontMANOVA::scMANOVA))

cat("semicontMANOVA version:", as.character(packageVersion("semicontMANOVA")), "\n")
cat("Installed at:", find.package("semicontMANOVA"), "\n")
cat("Contains missing.model argument:", has_missing_model, "\n")

if (!has_missing_model) {
  stop("The edited Bernoulli-enabled package was not installed correctly.")
}

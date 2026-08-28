###############################################################################
## Compare scMAN-MV and scMAN-Bernoulli on the Zeller application
###############################################################################

source(file.path("applications", "zeller_crc", "config.R"))
source(file.path("applications", "zeller_crc", "helpers.R"))
check_zeller_root()

prepared_path <- file.path(zeller_results_dir, "prepared_data.rds")
if (!file.exists(prepared_path)) stop("Run prepare_data.R first.")

if (!requireNamespace("semicontMANOVA", quietly = TRUE)) {
  stop("Install the edited thesis package with source('package/install_package.R').")
}
if (!("missing.model" %in% names(formals(semicontMANOVA::scMANOVA)))) {
  stop("The installed semicontMANOVA package is not the Bernoulli-enabled thesis version.")
}

prepared <- readRDS(prepared_path)
X <- prepared$X
group <- prepared$group
group_order <- prepared$group_order

row_order <- order(group)
X <- X[row_order, , drop = FALSE]
group <- droplevels(group[row_order])
n_vec <- as.numeric(table(group)[group_order])

num1 <- function(x) as.numeric(x[1])
retained_features <- function(fit) {
  if (!is.null(fit$mu)) return(ncol(fit$mu))
  if (!is.null(fit$mu0)) return(length(fit$mu0))
  if (!is.null(fit$sigma)) return(nrow(fit$sigma))
  NA_integer_
}
removed_features <- function(fit) {
  if (is.null(fit$removed.vars)) 0L else length(fit$removed.vars)
}

fit_one <- function(label, missing_model) {
  cat("Running", label, "with B =", B_perm, "...\n")
  set.seed(seed)
  start <- Sys.time()
  fit <- semicontMANOVA::scMANOVA(
    x = X,
    n = n_vec,
    lambda = lambda_bounds,
    lambda0 = lambda0_bounds,
    p.value.perm = p_value_perm,
    B = B_perm,
    ncpus = ncpus,
    parallel = parallel_mode,
    ident = ident,
    missing.model = missing_model
  )
  cat("Finished", label, "in",
      round(as.numeric(difftime(Sys.time(), start, units = "mins")), 1), "minutes.\n")

  ll_h1 <- num1(fit$logLik)
  ll_h0 <- num1(fit$logLik0)
  ll_miss_h1 <- num1(fit$logLikPi)
  ll_miss_h0 <- num1(fit$logLikPi0)
  ll_gauss_h1 <- ll_h1 - ll_miss_h1
  ll_gauss_h0 <- ll_h0 - ll_miss_h0

  stat_total <- num1(fit$statistic)
  stat_missing <- -2 * (ll_miss_h0 - ll_miss_h1)
  stat_gaussian <- -2 * (ll_gauss_h0 - ll_gauss_h1)

  result <- data.frame(
    model = label,
    missing_model = missing_model,
    n = nrow(X),
    p = ncol(X),
    group_sizes = paste(n_vec, collapse = ";"),
    B = B_perm,
    p_value = num1(fit$p.value),
    statistic = stat_total,
    lambda = num1(fit$lambda),
    lambda0 = num1(fit$lambda0),
    retained_gaussian_features = retained_features(fit),
    removed_gaussian_features = removed_features(fit)
  )

  decomposition <- data.frame(
    model = label,
    missing_model = missing_model,
    logLik_H1_total = ll_h1,
    logLik_H0_total = ll_h0,
    logLik_H1_missingness = ll_miss_h1,
    logLik_H0_missingness = ll_miss_h0,
    logLik_H1_gaussian = ll_gauss_h1,
    logLik_H0_gaussian = ll_gauss_h0,
    statistic_total = stat_total,
    statistic_missingness = stat_missing,
    statistic_gaussian = stat_gaussian,
    missingness_share_of_total = stat_missing / stat_total,
    gaussian_share_of_total = stat_gaussian / stat_total
  )

  list(result = result, decomposition = decomposition)
}

mv <- fit_one("scMAN-MV", "multivariate_bernoulli")
bern <- fit_one("scMAN-Bernoulli", "bernoulli")

write_csv(rbind(mv$result, bern$result), "model_results.csv")
write_csv(rbind(mv$decomposition, bern$decomposition), "likelihood_decomposition.csv")

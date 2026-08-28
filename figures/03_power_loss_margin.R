###############################################################################
## Thesis Figure: Gaussian-only versus full-test rejection margins
##
## Reproduces Figure \ref{fig:power-loss-margins} for the representative
## independent-Bernoulli mean-shift scenario:
##   c1 = 1, c2 = 0, pmiss = 0.5, gamma = 0, p = 50, n_k = 5.
##
## This is intentionally single-core and restartable. The compact intermediate
## CSV is written to figures/results/ and is ignored by Git.
###############################################################################

if (!requireNamespace("semicontMANOVA", quietly = TRUE)) {
  stop("Install the edited thesis package with source('package/install_package.R').")
}
if (!("missing.model" %in% names(formals(semicontMANOVA::scMANOVA)))) {
  stop("The installed semicontMANOVA package is not the Bernoulli-enabled thesis version.")
}

## Scenario and Monte Carlo settings used for the thesis diagnostic.
K <- 2L
n_k <- 5L
n_vec <- rep(n_k, K)
p <- 50L
c1 <- 1
pmiss <- 0.5
rho <- 0
nrep <- 250L
B <- 100L
alpha <- 0.05
seed_base <- 20260715L
lambda <- c(0, 100)
lambda0 <- c(0, 100)

models <- c("multivariate_bernoulli", "bernoulli")
model_labels <- c(
  multivariate_bernoulli = "Multivariate Bernoulli",
  bernoulli = "Independent Bernoulli"
)

result_dir <- file.path("figures", "results")
output_dir <- file.path("figures", "generated")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

result_file <- file.path(result_dir, "power_loss_margin_by_rep.csv")
figure_file <- file.path(output_dir, "power_loss_margin_scatter_by_model.png")

empty_results <- function() {
  data.frame(
    rep = integer(),
    model = character(),
    T_total = numeric(),
    T_gaussian = numeric(),
    p_total = numeric(),
    p_gaussian = numeric(),
    q95_total = numeric(),
    q95_gaussian = numeric(),
    full_margin = numeric(),
    gaussian_margin = numeric(),
    reject_full = logical(),
    reject_gaussian = logical(),
    lost_due_to_missingness = logical(),
    stringsAsFactors = FALSE
  )
}

read_results <- function() {
  if (!file.exists(result_file)) return(empty_results())
  x <- utils::read.csv(result_file, stringsAsFactors = FALSE)
  required <- names(empty_results())
  if (!all(required %in% names(x))) {
    stop("Unexpected columns in existing result file: ", result_file)
  }
  x[, required, drop = FALSE]
}

write_results <- function(x) {
  tmp <- paste0(result_file, ".tmp")
  utils::write.csv(x, tmp, row.names = FALSE)
  if (file.exists(result_file)) unlink(result_file)
  if (!file.rename(tmp, result_file)) {
    ok <- file.copy(tmp, result_file, overwrite = TRUE)
    unlink(tmp)
    if (!isTRUE(ok)) stop("Could not write: ", result_file)
  }
}

safe_num <- function(x) as.numeric(x[1])

fit_one_model <- function(x, missing_model) {
  semicontMANOVA::scMANOVA(
    x = x,
    n = n_vec,
    lambda = lambda,
    lambda0 = lambda0,
    p.value.perm = FALSE,
    B = 0L,
    ncpus = 1L,
    parallel = "no",
    ident = TRUE,
    missing.model = missing_model
  )
}

extract_components <- function(fit) {
  ll1 <- safe_num(fit$logLik)
  ll0 <- safe_num(fit$logLik0)
  ll1_miss <- safe_num(fit$logLikPi)
  ll0_miss <- safe_num(fit$logLikPi0)
  ll1_gauss <- ll1 - ll1_miss
  ll0_gauss <- ll0 - ll0_miss

  c(
    T_total = 2 * (ll1 - ll0),
    T_gaussian = 2 * (ll1_gauss - ll0_gauss)
  )
}

simulate_dataset <- function() {
  mu <- matrix(0, nrow = K, ncol = p)
  mu[2, ] <- c1
  semicontMANOVA::scMANOVAsimulation(
    n = n_vec,
    p = p,
    pmiss = rep(pmiss, K),
    rho = rho,
    mu = mu
  )
}

results <- read_results()

for (irep in seq_len(nrep)) {
  done_models <- results$model[results$rep == irep]
  missing_models <- setdiff(models, done_models)
  if (length(missing_models) == 0L) next

  ## This seed and the shared permutation indices reproduce the original
  ## thesis diagnostic. Both occurrence models see exactly the same data and
  ## exactly the same row permutations within each repetition.
  set.seed(seed_base + irep)
  x <- simulate_dataset()
  perm_indices <- replicate(B, sample.int(nrow(x)), simplify = FALSE)

  for (model_name in missing_models) {
    cat("repetition", irep, "/", nrep, "|", model_labels[[model_name]], "\n")
    fit_obs <- fit_one_model(x, model_name)
    obs <- extract_components(fit_obs)

    perm_total <- numeric(B)
    perm_gaussian <- numeric(B)

    for (b in seq_len(B)) {
      fit_perm <- fit_one_model(x[perm_indices[[b]], , drop = FALSE], model_name)
      comp <- extract_components(fit_perm)
      perm_total[b] <- comp[["T_total"]]
      perm_gaussian[b] <- comp[["T_gaussian"]]
    }

    p_total <- (1 + sum(perm_total >= obs[["T_total"]])) / (B + 1)
    p_gaussian <- (1 + sum(perm_gaussian >= obs[["T_gaussian"]])) / (B + 1)
    q95_total <- as.numeric(stats::quantile(perm_total, probs = 0.95, names = FALSE))
    q95_gaussian <- as.numeric(stats::quantile(perm_gaussian, probs = 0.95, names = FALSE))

    row <- data.frame(
      rep = irep,
      model = model_name,
      T_total = obs[["T_total"]],
      T_gaussian = obs[["T_gaussian"]],
      p_total = p_total,
      p_gaussian = p_gaussian,
      q95_total = q95_total,
      q95_gaussian = q95_gaussian,
      full_margin = obs[["T_total"]] - q95_total,
      gaussian_margin = obs[["T_gaussian"]] - q95_gaussian,
      reject_full = p_total <= alpha,
      reject_gaussian = p_gaussian <= alpha,
      lost_due_to_missingness = (p_gaussian <= alpha) && (p_total > alpha),
      stringsAsFactors = FALSE
    )

    results <- rbind(results, row)
    results <- results[order(results$rep, match(results$model, models)), , drop = FALSE]
    rownames(results) <- NULL
    write_results(results)
  }
}

if (nrow(results) != nrep * length(models)) {
  stop("The diagnostic is incomplete; expected ", nrep * length(models), " model/repetition rows.")
}

## Plot exactly the two model panels used in the thesis.
xlim_all <- range(c(results$gaussian_margin, 0), na.rm = TRUE)
ylim_all <- range(c(results$full_margin, 0), na.rm = TRUE)
x_pad <- 0.08 * diff(xlim_all)
y_pad <- 0.08 * diff(ylim_all)
if (!is.finite(x_pad) || x_pad == 0) x_pad <- 1
if (!is.finite(y_pad) || y_pad == 0) y_pad <- 1
xlim_all <- xlim_all + c(-x_pad, x_pad)
ylim_all <- ylim_all + c(-y_pad, y_pad)

plot_margin_panel <- function(d, main_title) {
  graphics::plot(
    d$gaussian_margin,
    d$full_margin,
    type = "n",
    xlim = xlim_all,
    ylim = ylim_all,
    xlab = "Gaussian-only margin",
    ylab = "Full-test margin",
    main = main_title,
    cex.lab = 1.55,
    cex.axis = 1.30,
    cex.main = 1.65
  )

  graphics::rect(
    xleft = 0,
    ybottom = ylim_all[1],
    xright = xlim_all[2],
    ytop = 0,
    col = grDevices::rgb(1, 0, 0, alpha = 0.08),
    border = NA
  )

  graphics::grid(col = "grey90")
  graphics::abline(v = 0, lty = 2, lwd = 1.3, col = "grey35")
  graphics::abline(h = 0, lty = 2, lwd = 1.3, col = "grey35")

  point_col <- ifelse(
    d$lost_due_to_missingness,
    "#D55E00",
    ifelse(d$reject_full, "#0072B2", "#666666")
  )
  point_pch <- ifelse(d$lost_due_to_missingness, 19, ifelse(d$reject_full, 17, 1))

  graphics::points(
    d$gaussian_margin,
    d$full_margin,
    pch = point_pch,
    col = point_col,
    bg = point_col,
    cex = 1.10
  )

  graphics::legend(
    "topleft",
    legend = c(
      "Lost: Gaussian rejects, full does not",
      "Full rejects",
      "No full rejection"
    ),
    pch = c(19, 17, 1),
    col = c("#D55E00", "#0072B2", "#666666"),
    bty = "n",
    cex = 1.05
  )

  graphics::text(
    x = xlim_all[2],
    y = ylim_all[1],
    labels = "Gaussian signal present\nbut full test fails",
    adj = c(1.02, -0.15),
    cex = 1.05,
    col = "#D55E00"
  )
}

grDevices::png(
  filename = figure_file,
  width = 3200,
  height = 1200,
  res = 180
)
graphics::par(mfrow = c(1, 2), mar = c(5.4, 5.4, 4.4, 1.2) + 0.1)
for (model_name in models) {
  d <- results[results$model == model_name, , drop = FALSE]
  plot_margin_panel(d, model_labels[[model_name]])
}
graphics::par(mfrow = c(1, 1))
grDevices::dev.off()

cat("Wrote thesis figure:", figure_file, "\n")

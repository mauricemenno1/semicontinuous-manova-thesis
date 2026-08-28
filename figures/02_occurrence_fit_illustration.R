###############################################################################
## Single-dataset illustration of sparsity in the occurrence models
##
## Setting matched to the multivariate-missingness simulation study:
##   p   = 50 components
##   K   = 2 groups
##   n_k = 5 observations per group
##   H1_2 probability shift: pmiss_1 = 0.50, pmiss_2 = 0.65
##   exchangeable-indicator ICC = 0.05
##
## Two occurrence mechanisms are generated:
##   1. Independent Bernoulli indicators for every component.
##   2. Exchangeable multivariate-Bernoulli indicators generated count-first:
##      S follows the same beta-binomial distribution as in
##      simulate_missingness_mv/01_simulation_helpers.R; conditional on S = s,
##      s positions are selected uniformly without replacement.
##
## Plotting:
##   - Bernoulli feature-specific estimates are shown as points only.
##   - Bernoulli-implied count distributions are shown as continuous curves.
##   - Exchangeable count estimates are shown as vertical spikes at the
##     observed count values.
##
## Exactly four thesis PNG files are produced (the group-1 panels used in the thesis).
###############################################################################

## User-adjustable settings ----------------------------------------------------

seed <- 20260728L

K <- 2L
n_k <- 5L
n_vec <- rep(n_k, K)
p <- 50L

## H1_2 probability-shift setting used in the simulation study.
base_pmiss <- 0.50
c2 <- 0.15
pmiss_by_group <- base_pmiss + c2 * seq(0, K - 1) / (K - 1)
presence_probability <- 1 - pmiss_by_group

## Exchangeable occurrence correlation.
mv_missingness_icc <- 0.05

## Output folder.
output_dir <- file.path("figures", "generated", "occurrence_fit_illustration")

## High-resolution thesis-ready PNG dimensions.
plot_width_px <- 2000
plot_height_px <- 1650L
plot_resolution <- 300L

## Colours.
true_colour <- "blue"
estimate_colour <- "red"

## Maximum y-axis value for all count-distribution plots.
common_count_y_max <- 0.20

## Checks ----------------------------------------------------------------------

if (K != 2L) stop("This illustration is written for K = 2 groups.")
if (length(n_vec) != K) stop("n_vec must contain one sample size per group.")

if (length(pmiss_by_group) != K) {
  stop("pmiss_by_group must contain one value per group.")
}

if (any(pmiss_by_group <= 0 | pmiss_by_group >= 1)) {
  stop("All missingness probabilities must lie strictly between 0 and 1.")
}

if (!is.finite(mv_missingness_icc) ||
    mv_missingness_icc <= 0 ||
    mv_missingness_icc >= 1) {
  stop("mv_missingness_icc must lie strictly between 0 and 1.")
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

## Remove old PNGs from this dedicated output folder.
old_pngs <- list.files(
  output_dir,
  pattern = "\\.png$",
  full.names = TRUE
)

if (length(old_pngs) > 0L) {
  unlink(old_pngs)
}

## Helper functions ------------------------------------------------------------

## Beta-binomial count distribution.
positive_count_distribution <- function(p, pmiss, icc) {
  
  p <- as.integer(p)
  
  probabilities <- numeric(p + 1L)
  names(probabilities) <- as.character(0:p)
  
  if (pmiss <= 0) {
    probabilities[p + 1L] <- 1
    return(probabilities)
  }
  
  if (pmiss >= 1) {
    probabilities[1L] <- 1
    return(probabilities)
  }
  
  concentration <- (1 / icc) - 1
  
  shape_positive <- (1 - pmiss) * concentration
  shape_zero <- pmiss * concentration
  
  s <- 0:p
  
  log_probabilities <-
    lchoose(p, s) +
    lbeta(
      s + shape_positive,
      p - s + shape_zero
    ) -
    lbeta(
      shape_positive,
      shape_zero
    )
  
  log_probabilities <-
    log_probabilities - max(log_probabilities)
  
  probabilities <- exp(log_probabilities)
  probabilities <- probabilities / sum(probabilities)
  
  names(probabilities) <- as.character(s)
  
  probabilities
}

## Independent component-wise Bernoulli occurrence generation.
simulate_independent_presence_patterns <- function(
    n_obs,
    p,
    pmiss
) {
  
  matrix(
    stats::rbinom(
      n = n_obs * p,
      size = 1L,
      prob = 1 - pmiss
    ),
    nrow = n_obs,
    ncol = p
  )
}

## Count-first exchangeable generator.
simulate_exchangeable_presence_patterns <- function(
    n_obs,
    p,
    pmiss,
    icc
) {
  
  count_probabilities <- positive_count_distribution(
    p = p,
    pmiss = pmiss,
    icc = icc
  )
  
  ## Step 1: draw S, the number of positive components.
  positive_counts <- sample.int(
    n = p + 1L,
    size = n_obs,
    replace = TRUE,
    prob = count_probabilities
  ) - 1L
  
  ## Step 2: conditional on S = s, choose s positions uniformly.
  y <- matrix(
    0L,
    nrow = n_obs,
    ncol = p
  )
  
  for (i in seq_len(n_obs)) {
    
    s <- positive_counts[i]
    
    if (s > 0L) {
      
      positive_positions <- sample.int(
        p,
        size = s,
        replace = FALSE
      )
      
      y[i, positive_positions] <- 1L
    }
  }
  
  y
}

## Distribution of the sum of independent, non-identically distributed
## Bernoulli variables.
poisson_binomial_pmf <- function(probabilities) {
  
  probabilities <- as.numeric(probabilities)
  m <- length(probabilities)
  
  probability <- numeric(m + 1L)
  probability[1L] <- 1
  
  for (rho in probabilities) {
    
    old <- probability
    
    probability <- old * (1 - rho)
    
    probability[2:(m + 1L)] <-
      probability[2:(m + 1L)] +
      old[1:m] * rho
  }
  
  probability[
    probability < 0 &
      probability > -1e-14
  ] <- 0
  
  probability / sum(probability)
}

## Group-specific exchangeable count-model MLE:
## q_k(s) = P(S = s).
empirical_count_pmf <- function(y, p) {
  
  tabulate(
    rowSums(y) + 1L,
    nbins = p + 1L
  ) / nrow(y)
}

## Open high-resolution PNG.
open_png <- function(filename) {
  
  grDevices::png(
    filename = file.path(
      output_dir,
      filename
    ),
    width = plot_width_px,
    height = plot_height_px,
    res = plot_resolution
  )
  
  graphics::par(
    mar = c(5.2, 5.2, 1.2, 1.2),
    mgp = c(3.2, 0.9, 0),
    las = 1,
    cex.axis = 1.15,
    cex.lab = 1.25,
    family = "serif",
    bty = "o",
    xaxs = "r",
    yaxs = "r"
  )
}

close_png <- function() {
  grDevices::dev.off()
}

## Legend for the feature-specific Bernoulli plots.
add_feature_legend <- function() {
  
  graphics::legend(
    "topleft",
    legend = c(
      "True presence probability",
      "Bernoulli estimate"
    ),
    col = c(
      true_colour,
      estimate_colour
    ),
    lty = c(1, NA),
    lwd = c(2.5, NA),
    pch = c(NA, 16),
    pt.cex = 0.80,
    cex = 1.05,
    bty = "n"
  )
}

## Legend for count-distribution plots.
add_count_legend <- function(
    estimate_label,
    estimate_style
) {
  
  if (estimate_style == "curve") {
    
    graphics::legend(
      "topleft",
      legend = c(
        "True distribution of S",
        estimate_label
      ),
      col = c(
        true_colour,
        estimate_colour
      ),
      lty = c(1, 1),
      lwd = c(2.5, 2.2),
      pch = c(NA, NA),
      cex = 1.05,
      bty = "n"
    )
    
  } else if (estimate_style == "spikes") {
    
    graphics::legend(
      "topleft",
      legend = c(
        "True distribution of S",
        estimate_label
      ),
      col = c(
        true_colour,
        estimate_colour
      ),
      lty = c(1, NA),
      lwd = c(2.5, NA),
      pch = c(NA, 124),
      pt.cex = 1.3,
      cex = 1.05,
      bty = "n"
    )
  }
}

## Feature-specific Bernoulli probabilities.
##
## The true probability remains a blue horizontal line.
## The fitted probabilities are red points only.
plot_feature_probabilities <- function(
    estimates,
    true_probability,
    filename
) {
  
  component <- seq_along(estimates)
  
  open_png(filename)
  
  graphics::plot(
    component,
    estimates,
    type = "n",
    xlim = c(
      1,
      length(estimates)
    ),
    ylim = c(0, 1),
    xlab = "Component index, j",
    ylab = "Presence probability"
  )
  
  graphics::abline(
    h = true_probability,
    col = true_colour,
    lwd = 2.5
  )
  
  graphics::points(
    component,
    estimates,
    col = estimate_colour,
    pch = 16,
    cex = 0.70
  )
  
  add_feature_legend()
  
  close_png()
}

## Plot a distribution of S.
##
## estimate_style = "curve":
##   fitted distribution is drawn as a connected red curve.
##
## estimate_style = "spikes":
##   only the nonzero fitted probabilities are shown, with a straight
##   vertical line from zero to each probability mass.
plot_count_distribution <- function(
    true_distribution,
    estimated_distribution,
    filename,
    estimate_label,
    estimate_style,
    y_max
) {
  
  s <- 0:(length(true_distribution) - 1L)
  
  open_png(filename)
  
  graphics::plot(
    s,
    true_distribution,
    type = "n",
    xlim = range(s),
    ylim = c(0, y_max),
    xlab = "Number of positive components, s",
    ylab = "Probability"
  )
  
  ## True distribution remains unchanged.
  graphics::lines(
    s,
    true_distribution,
    col = true_colour,
    lwd = 2.5
  )
  
  if (estimate_style == "curve") {
    
    ## Bernoulli-implied distribution:
    ## retain the connected red curve.
    graphics::lines(
      s,
      estimated_distribution,
      col = estimate_colour,
      lwd = 2.2
    )
    
  } else if (estimate_style == "spikes") {
    
    ## Exchangeable count estimate:
    ## straight vertical lines at the observed count values.
    nonzero <- estimated_distribution > 0
    
    graphics::segments(
      x0 = s[nonzero],
      y0 = 0,
      x1 = s[nonzero],
      y1 = estimated_distribution[nonzero],
      col = estimate_colour,
      lwd = 2.2
    )
    
    ## Retain a small point at the top of each probability mass.
    graphics::points(
      s[nonzero],
      estimated_distribution[nonzero],
      col = estimate_colour,
      pch = 16,
      cex = 0.80
    )
  }
  
  add_count_legend(
    estimate_label = estimate_label,
    estimate_style = estimate_style
  )
  
  close_png()
}

## Generate one group-specific dataset under each mechanism -------------------

set.seed(seed)

## Independent Bernoulli injection.
y_independent_by_group <- lapply(
  seq_len(K),
  function(k) {
    
    simulate_independent_presence_patterns(
      n_obs = n_vec[k],
      p = p,
      pmiss = pmiss_by_group[k]
    )
  }
)

## Exchangeable beta-binomial injection.
y_multivariate_by_group <- lapply(
  seq_len(K),
  function(k) {
    
    simulate_exchangeable_presence_patterns(
      n_obs = n_vec[k],
      p = p,
      pmiss = pmiss_by_group[k],
      icc = mv_missingness_icc
    )
  }
)

## True group-specific distributions of S -------------------------------------

s_values <- 0:p

## Independent Bernoulli generator:
## S follows a binomial distribution.
true_q_independent_by_group <- lapply(
  seq_len(K),
  function(k) {
    
    stats::dbinom(
      s_values,
      size = p,
      prob = presence_probability[k]
    )
  }
)

## Exchangeable beta-binomial generator.
true_q_multivariate_by_group <- lapply(
  seq_len(K),
  function(k) {
    
    positive_count_distribution(
      p = p,
      pmiss = pmiss_by_group[k],
      icc = mv_missingness_icc
    )
  }
)

## Pre-compute all group-specific H1 estimates --------------------------------

## Independent injection fitted with feature-specific Bernoulli model.
rho_hat_independent_by_group <- lapply(
  seq_len(K),
  function(k) {
    colMeans(y_independent_by_group[[k]])
  }
)

## Independent injection fitted with exchangeable count model.
q_hat_independent_mv_by_group <- lapply(
  seq_len(K),
  function(k) {
    
    empirical_count_pmf(
      y = y_independent_by_group[[k]],
      p = p
    )
  }
)

## Exchangeable injection fitted with Bernoulli model.
rho_hat_multivariate_by_group <- lapply(
  seq_len(K),
  function(k) {
    colMeans(y_multivariate_by_group[[k]])
  }
)

## Count distribution implied by fitted Bernoulli probabilities.
q_hat_bernoulli_implied_by_group <- lapply(
  seq_len(K),
  function(k) {
    
    poisson_binomial_pmf(
      rho_hat_multivariate_by_group[[k]]
    )
  }
)

## Exchangeable injection fitted with exchangeable count model.
q_hat_multivariate_mv_by_group <- lapply(
  seq_len(K),
  function(k) {
    
    empirical_count_pmf(
      y = y_multivariate_by_group[[k]],
      p = p
    )
  }
)

## Group-1 H1 fits and four thesis plots ---------------------------------------
##
## Both groups are still simulated above. This deliberately preserves the RNG
## sequence of the original plotting script, so the retained group-1 exchangeable
## panels are identical to the figures used in the thesis. Only group 1 is plotted.

for (k in 1L) {
  
  group_tag <- paste0(
    "h1_group_",
    k
  )
  
  ## 1. Independent injection fitted by feature-specific Bernoulli model.
  ##
  ## Red estimates are POINTS ONLY.
  rho_hat_independent <- rho_hat_independent_by_group[[k]]
  
  plot_feature_probabilities(
    estimates = rho_hat_independent,
    true_probability = presence_probability[k],
    filename = paste0(
      "independent_injection_bernoulli_fit_",
      group_tag,
      "_feature_presence_probabilities.png"
    )
  )
  
  ## 2. Independent injection fitted by exchangeable count model.
  ##
  ## Red estimate is shown as vertical probability-mass SPIKES.
  q_hat_independent_mv <- q_hat_independent_mv_by_group[[k]]
  
  plot_count_distribution(
    true_distribution = true_q_independent_by_group[[k]],
    estimated_distribution = q_hat_independent_mv,
    filename = paste0(
      "independent_injection_multivariate_fit_",
      group_tag,
      "_count_distribution.png"
    ),
    estimate_label = "Multivariate count estimate",
    estimate_style = "spikes",
    y_max = common_count_y_max
  )
  
  ## 3. Exchangeable beta-binomial injection fitted by Bernoulli model.
  ##
  ## Bernoulli-implied count distribution remains a connected red CURVE.
  q_hat_bernoulli_implied <-
    q_hat_bernoulli_implied_by_group[[k]]
  
  plot_count_distribution(
    true_distribution = true_q_multivariate_by_group[[k]],
    estimated_distribution = q_hat_bernoulli_implied,
    filename = paste0(
      "beta_binomial_injection_bernoulli_fit_",
      group_tag,
      "_implied_count_distribution.png"
    ),
    estimate_label = "Bernoulli-implied distribution",
    estimate_style = "curve",
    y_max = common_count_y_max
  )
  
  ## 4. Exchangeable beta-binomial injection fitted by exchangeable count model.
  ##
  ## Red estimate is shown as vertical probability-mass SPIKES.
  q_hat_multivariate_mv <- q_hat_multivariate_mv_by_group[[k]]
  
  plot_count_distribution(
    true_distribution = true_q_multivariate_by_group[[k]],
    estimated_distribution = q_hat_multivariate_mv,
    filename = paste0(
      "beta_binomial_injection_multivariate_fit_",
      group_tag,
      "_count_distribution.png"
    ),
    estimate_label = "Multivariate count estimate",
    estimate_style = "spikes",
    y_max = common_count_y_max
  )
}

## Final check -----------------------------------------------------------------

created_pngs <- list.files(
  output_dir,
  pattern = "\\.png$",
  full.names = TRUE
)

if (length(created_pngs) != 4L) {
  
  stop(
    "Expected exactly four PNG files, but found ",
    length(created_pngs),
    "."
  )
}

cat(
  "Finished. Four group-1 thesis PNG figures were written to: ",
  normalizePath(output_dir),
  "\n",
  sep = ""
)
###############################################################################
## Thesis Figure: binomial versus beta-binomial count distribution
##
## Reproduces Figure \ref{fig:distribution-s} from the thesis.
## Run from the repository root.
###############################################################################

output_dir <- file.path("figures", "generated")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

p <- 50L
zero_probability <- 0.5
dependence <- 0.05
positive_probability <- 1 - zero_probability
s <- 0:p

prob_binomial <- stats::dbinom(s, size = p, prob = positive_probability)

concentration <- (1 / dependence) - 1
alpha <- positive_probability * concentration
beta <- zero_probability * concentration

log_prob_beta_binomial <-
  lchoose(p, s) +
  lbeta(s + alpha, p - s + beta) -
  lbeta(alpha, beta)

prob_beta_binomial <- exp(log_prob_beta_binomial)
prob_beta_binomial <- prob_beta_binomial / sum(prob_beta_binomial)

draw_distribution_plot <- function() {
  old_par <- graphics::par(
    mar = c(5.2, 5.2, 1.2, 1.2),
    mgp = c(3.2, 0.9, 0),
    las = 1,
    cex.axis = 1.15,
    cex.lab = 1.25
  )
  on.exit(graphics::par(old_par))

  graphics::plot(
    s,
    prob_binomial,
    type = "l",
    lty = 1,
    lwd = 2.5,
    col = "blue",
    xlab = "Number of positive components, s",
    ylab = "Probability",
    ylim = c(0, 1.05 * max(prob_binomial, prob_beta_binomial))
  )

  graphics::lines(
    s,
    prob_beta_binomial,
    lty = 1,
    lwd = 2.5,
    col = "red"
  )

  graphics::legend(
    "topright",
    legend = c(
      "Independent Bernoulli: binomial",
      "Exchangeable model: beta-binomial"
    ),
    col = c("blue", "red"),
    lty = c(1, 1),
    lwd = c(2.5, 2.5),
    cex = 1.05,
    bty = "n"
  )
}

png_args <- list(
  filename = file.path(output_dir, "distribution_S_binomial_vs_beta_binomial.png"),
  width = 2800,
  height = 1650,
  res = 300
)

## The thesis PNG was generated with Cairo. Standard Windows and Linux R builds
## normally support it; the fallback still produces the same plotted quantities.
if (isTRUE(capabilities("cairo"))) png_args$type <- "cairo"
do.call(grDevices::png, png_args)
draw_distribution_plot()
grDevices::dev.off()

cat("Wrote:", png_args$filename, "\n")

## -----------------------------------------------------------------------------
## Custom missingness plots for the Zeller dataset
##
## Run from the repository root with:
##   source("applications/zeller_crc/make_figures.R")
##
## Output folder:
##   applications/zeller_crc/figures/
##
## This script generates only the four selected thesis figures.
## -----------------------------------------------------------------------------

## ==============================
## User-adjustable settings
## ==============================

# Number of largest absolute group differences shown in the Bernoulli plots.
top_n_features <- 50

# Bin width for the count-based model.
count_bin_width <- 5

# Colours.
col_ctr <- "#2C7BB6"
col_crc <- "#D7191C"
col_h0  <- "#4D4D4D"

# Figure size and resolution.
png_width  <- 1800
png_height <- 1200
png_res    <- 220

# Easily adjustable text sizes for LaTeX use.
axis_text_cex   <- 1.55
axis_label_cex  <- 1.75
legend_cex      <- 1.45
line_width_main <- 3.6
line_width_aux  <- 2.6
point_cex       <- 1.2

# Margins.
plot_mar <- c(5.4, 6.0, 1.4, 2.0)

# Small expansion factor around data-driven y-ranges.
ylim_expand <- 0.08

## ==============================
## Helper functions
## ==============================

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
msg <- function(...) cat(timestamp(), "|", ..., "\n")

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_\\.-]", "_", x)
}

range_expand <- function(x, frac = 0.08, include_zero = FALSE) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(0, 1))
  if (include_zero) x <- c(x, 0)
  r <- range(x)
  if (r[1] == r[2]) {
    bump <- ifelse(r[1] == 0, 0.01, abs(r[1]) * frac)
    return(c(r[1] - bump, r[2] + bump))
  }
  d <- diff(r)
  c(r[1] - frac * d, r[2] + frac * d)
}

make_bin_prob <- function(counts, breaks, bin_levels) {
  bins <- cut(counts, breaks = breaks, include.lowest = TRUE, right = FALSE)
  tab <- table(factor(bins, levels = bin_levels))
  as.numeric(tab) / length(counts)
}

set_plot_par <- function() {
  par(
    mar = plot_mar,
    cex.axis = axis_text_cex,
    cex.lab = axis_label_cex,
    family = "sans"
  )
}

open_png <- function(filename) {
  png(filename = filename, width = png_width, height = png_height, res = png_res)
}

## ==============================
## Locate files
## ==============================

base_dir <- file.path("applications", "zeller_crc")
out_dir <- file.path(base_dir, "figures")
input_path <- file.path(base_dir, "results", "prepared_data.rds")

if (!file.exists(file.path(base_dir, "config.R"))) {
  stop("Run this script from the repository root.")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) {
  stop(
    "Could not find prepared data at:\n", input_path,
    "\nFirst run source('applications/zeller_crc/prepare_data.R')."
  )
}

## Remove old PNGs so this folder contains only the four thesis figures.
old_pngs <- list.files(out_dir, pattern = "\\.png$", full.names = TRUE)
if (length(old_pngs) > 0L) unlink(old_pngs)

## ==============================
## Load prepared data
## ==============================

msg("Loading prepared Zeller data...")
prepared <- readRDS(input_path)

X <- prepared$X
group <- droplevels(prepared$group)

if (!is.null(prepared$group_order)) {
  group <- factor(as.character(group), levels = prepared$group_order)
}

groups <- levels(group)
if (length(groups) != 2L) {
  stop("Expected exactly two groups, found: ", paste(groups, collapse = ", "))
}

g1 <- groups[1]
g2 <- groups[2]

Y <- X > 0

msg("Loaded data with ", nrow(X), " samples and ", ncol(X), " features.")
msg("Group counts: ", paste(names(table(group)), as.integer(table(group)), collapse = "; "))

## ==============================
## Count-based quantities
## ==============================

positive_count <- rowSums(Y)
S1 <- positive_count[group == g1]
S2 <- positive_count[group == g2]
S0 <- positive_count

break_min <- floor(min(positive_count) / count_bin_width) * count_bin_width
break_max <- ceiling((max(positive_count) + 1) / count_bin_width) * count_bin_width
breaks <- seq(break_min, break_max, by = count_bin_width)
bin_mid <- head(breaks, -1) + diff(breaks) / 2
bin_levels <- levels(cut(positive_count, breaks = breaks, include.lowest = TRUE, right = FALSE))

bpmf_g1 <- make_bin_prob(S1, breaks, bin_levels)
bpmf_g2 <- make_bin_prob(S2, breaks, bin_levels)
bpmf_h0 <- make_bin_prob(S0, breaks, bin_levels)

diff_count_g1_h0 <- bpmf_g1 - bpmf_h0
diff_count_g2_h0 <- bpmf_g2 - bpmf_h0

count_ylim <- range_expand(c(bpmf_g1, bpmf_g2, bpmf_h0), frac = ylim_expand, include_zero = TRUE)
count_diff_ylim <- range_expand(c(diff_count_g1_h0, diff_count_g2_h0), frac = 0.10, include_zero = TRUE)

## ==============================
## Bernoulli feature-wise quantities
## ==============================

rho_g1 <- colMeans(Y[group == g1, , drop = FALSE])
rho_g2 <- colMeans(Y[group == g2, , drop = FALSE])
rho_h0 <- colMeans(Y)

signed_diff <- rho_g2 - rho_g1
abs_diff <- abs(signed_diff)
# Final thesis ordering: largest absolute CTR-CRC prevalence difference first.
feature_order <- order(abs_diff, decreasing = TRUE)

rho1_ord <- rho_g1[feature_order]
rho2_ord <- rho_g2[feature_order]
rho0_ord <- rho_h0[feature_order]

# Differences relative to H0.
diff_rho1_h0 <- rho1_ord - rho0_ord
diff_rho2_h0 <- rho2_ord - rho0_ord

x_idx <- seq_along(rho1_ord)

if (length(top_n_features) != 1L || !is.finite(top_n_features) || top_n_features < 1) {
  stop("top_n_features must be one positive number.")
}

top_n_use <- min(as.integer(top_n_features), length(x_idx))
top_idx <- seq_len(top_n_use)

# Top features by absolute group difference.
rho1_top_raw <- rho1_ord[top_idx]
rho2_top_raw <- rho2_ord[top_idx]
rho0_top_raw <- rho0_ord[top_idx]

top_raw_ylim <- range_expand(c(rho1_top_raw, rho2_top_raw, rho0_top_raw), frac = 0.06)
top_raw_ylim[1] <- max(0, top_raw_ylim[1])
top_raw_ylim[2] <- min(1, top_raw_ylim[2])

# For the H1-H0 difference plot, use ALL features. This plot is deliberately
# independent of top_n_features: every feature in the dataset is included.
# Reorder by the CRC difference from H0: positive differences from largest to
# smallest, followed by negative differences from closest to zero to most
# negative. The CTR estimates are then matched to the same feature order.
all_feature_idx <- seq_along(diff_rho2_h0)

positive_features <- all_feature_idx[diff_rho2_h0 > 0]
zero_features     <- all_feature_idx[diff_rho2_h0 == 0]
negative_features <- all_feature_idx[diff_rho2_h0 < 0]

positive_features <- positive_features[
  order(diff_rho2_h0[positive_features], decreasing = TRUE)
]
negative_features <- negative_features[
  order(diff_rho2_h0[negative_features], decreasing = TRUE)
]

all_reordered_idx <- c(positive_features, zero_features, negative_features)

if (length(all_reordered_idx) != length(diff_rho2_h0)) {
  stop("The all-feature Bernoulli difference plot does not contain every feature.")
}

n_all_features <- length(all_reordered_idx)
if (n_all_features != ncol(Y)) {
  stop("Expected the all-feature plot to contain ", ncol(Y),
       " features, but it contains ", n_all_features, ".")
}

plot6_x <- seq_along(all_reordered_idx)
plot6_crc <- diff_rho2_h0[all_reordered_idx]
plot6_ctr <- diff_rho1_h0[all_reordered_idx]
plot6_ylim <- range_expand(c(plot6_crc, plot6_ctr), frac = 0.08, include_zero = TRUE)


## ==============================
## Plot 1
## Original count-based model: H1 and H0 fitted count distributions
## ==============================

file1 <- file.path(out_dir, safe_filename("Original_count_based_model_H1_and_H0_fitted_count_distributions.png"))
open_png(file1)
set_plot_par()
plot(
  bin_mid, bpmf_h0,
  type = "b", pch = 15, lwd = line_width_main, cex = point_cex, col = col_h0,
  xlab = paste0("Number of positive features per sample (bin width = ", count_bin_width, ")"),
  ylab = "Fitted probability",
  ylim = count_ylim
)
lines(bin_mid, bpmf_g1, type = "b", pch = 16, lwd = line_width_main, cex = point_cex, col = col_ctr)
lines(bin_mid, bpmf_g2, type = "b", pch = 17, lwd = line_width_main, cex = point_cex, col = col_crc)
grid()
legend(
  "topright",
  legend = c("H0 pooled", paste0("H1 ", g1), paste0("H1 ", g2)),
  col = c(col_h0, col_ctr, col_crc),
  pch = c(15, 16, 17),
  lwd = c(line_width_main, line_width_main, line_width_main),
  pt.cex = point_cex,
  cex = legend_cex,
  bty = "n"
)
dev.off()
msg("Wrote: ", file1)

## ==============================
## Plot 2
## Original count-based model: group PMFs minus H0 pooled PMF
## ==============================

file2 <- file.path(out_dir, safe_filename("Original_count_based_model_group_PMFs_minus_H0_pooled_PMF.png"))
open_png(file2)
set_plot_par()
plot(
  bin_mid, diff_count_g1_h0,
  type = "b", pch = 16, lwd = line_width_aux, cex = point_cex, col = col_ctr,
  xlab = paste0("Number of positive features per sample (bin width = ", count_bin_width, ")"),
  ylab = "Fitted probability difference",
  ylim = count_diff_ylim
)
abline(h = 0, lty = 2, lwd = 1.8, col = "grey40")
lines(bin_mid, diff_count_g2_h0, type = "b", pch = 17, lwd = line_width_aux, cex = point_cex, col = col_crc)
grid()
legend(
  "topright",
  legend = c(paste0("H1 ", g1, " - H0"), paste0("H1 ", g2, " - H0")),
  col = c(col_ctr, col_crc),
  pch = c(16, 17),
  lwd = c(line_width_aux, line_width_aux),
  pt.cex = point_cex,
  cex = legend_cex,
  bty = "n"
)
dev.off()
msg("Wrote: ", file2)

## ==============================
## Plot 3
## Bernoulli model under H0 and H1:
## raw fitted probabilities for the top features.
## ==============================

file3 <- file.path(
  out_dir,
  safe_filename(
    paste0(
      "Bernoulli_model_top_",
      top_n_use,
      "_raw_feature_wise_fitted_probabilities_H0_CTR_CRC.png"
    )
  )
)

open_png(file3)
set_plot_par()
plot(
  top_idx, rho0_top_raw,
  type = "n",
  xlim = c(0.5, top_n_use + 0.5),
  ylim = top_raw_ylim,
  xlab = paste0("Top ", top_n_use, " features, ordered by absolute fitted\n", "presence-probability difference (CTR vs CRC)"),
  ylab = "Fitted presence probability",
  xaxt = "n"
)
segments(
  x0 = top_idx,
  y0 = pmin(rho1_top_raw, rho2_top_raw, rho0_top_raw),
  x1 = top_idx,
  y1 = pmax(rho1_top_raw, rho2_top_raw, rho0_top_raw),
  col = adjustcolor("grey55", 0.65),
  lwd = 1.5
)
points(top_idx - 0.18, rho1_top_raw, pch = 16, cex = point_cex, col = col_ctr)
points(top_idx, rho0_top_raw, pch = 15, cex = point_cex, col = col_h0)
points(top_idx + 0.18, rho2_top_raw, pch = 17, cex = point_cex, col = col_crc)
grid()
legend(
  "topright",
  legend = c(paste0("H1 ", g1), "H0 pooled", paste0("H1 ", g2)),
  col = c(col_ctr, col_h0, col_crc),
  pch = c(16, 15, 17),
  pt.cex = point_cex,
  cex = legend_cex,
  bty = "n"
)
dev.off()
msg("Wrote: ", file3)

## ==============================
## Plot 4
## Bernoulli model: group probabilities minus H0 pooled probabilities
## using all features, reordered by the CRC difference from H0.
## Dots only, no connecting lines.
## ==============================

msg("Bernoulli H1-H0 plot uses ALL ", n_all_features, " features. top_n_features = ", top_n_features, " is ignored for this plot.")
file4 <- file.path(out_dir, safe_filename("Bernoulli_model_group_probabilities_minus_H0_pooled_probabilities.png"))
open_png(file4)
set_plot_par()
plot(
  plot6_x, plot6_crc,
  type = "n",
  xlim = c(0.5, n_all_features + 0.5),
  ylim = plot6_ylim,
  xlab = paste0("All ", n_all_features, " features, ordered by fitted\n",
                "CRC - H0 presence-probability difference"),
  ylab = "Fitted probability difference from H0",
  xaxt = "n"
)
abline(h = 0, lty = 2, lwd = 1.8, col = "grey40")
points(plot6_x, plot6_ctr, pch = 16, cex = 0.75 * point_cex, col = col_ctr)
points(plot6_x, plot6_crc, pch = 17, cex = 0.75 * point_cex, col = col_crc)
grid()
legend(
  "topright",
  legend = c(paste0(g1, " - H0"), paste0(g2, " - H0")),
  col = c(col_ctr, col_crc),
  pch = c(16, 17),
  pt.cex = point_cex,
  cex = legend_cex,
  bty = "n"
)
dev.off()
msg("Wrote ALL-FEATURE plot with ", n_all_features, " features: ", file4)

created_pngs <- list.files(out_dir, pattern = "\\.png$", full.names = TRUE)
if (length(created_pngs) != 4L) {
  stop("Expected exactly four PNG files, but found ", length(created_pngs), ".")
}

msg("Done. Selected images written to: ", out_dir)

###############################################################################
## Descriptive summaries used in the Zeller thesis application
###############################################################################

source(file.path("applications", "zeller_crc", "config.R"))
source(file.path("applications", "zeller_crc", "helpers.R"))
check_zeller_root()

prepared_path <- file.path(zeller_results_dir, "prepared_data.rds")
if (!file.exists(prepared_path)) stop("Run prepare_data.R first.")

prepared <- readRDS(prepared_path)
X <- prepared$X
group <- prepared$group

feature_zero <- colMeans(X == 0)
sample_zero <- rowMeans(X == 0)
sample_positive_count <- rowSums(X > 0)

q <- function(x, prob) as.numeric(stats::quantile(x, probs = prob, names = FALSE))

write_csv(
  data.frame(
    metric = c(
      "n_samples", "n_features", "overall_zero_fraction", "overall_positive_fraction",
      "feature_zero_min", "feature_zero_q25", "feature_zero_median", "feature_zero_mean",
      "feature_zero_q75", "feature_zero_max", "feature_zero_iqr", "features_zero_gt_50_pct",
      "features_zero_gt_80_pct", "sample_zero_min", "sample_zero_q25", "sample_zero_median",
      "sample_zero_mean", "sample_zero_q75", "sample_zero_max", "sample_positive_count_median"
    ),
    value = c(
      nrow(X), ncol(X), mean(X == 0), mean(X > 0),
      min(feature_zero), q(feature_zero, 0.25), median(feature_zero), mean(feature_zero),
      q(feature_zero, 0.75), max(feature_zero), stats::IQR(feature_zero), mean(feature_zero > 0.50),
      mean(feature_zero > 0.80), min(sample_zero), q(sample_zero, 0.25), median(sample_zero),
      mean(sample_zero), q(sample_zero, 0.75), max(sample_zero), median(sample_positive_count)
    )
  ),
  "dataset_sparsity_summary.csv"
)

group_summary <- do.call(rbind, lapply(group_order, function(g) {
  Xg <- X[group == g, , drop = FALSE]
  data.frame(
    group = g,
    n = nrow(Xg),
    p = ncol(Xg),
    overall_zero_fraction = mean(Xg == 0),
    overall_positive_fraction = mean(Xg > 0),
    median_feature_zero_fraction = median(colMeans(Xg == 0)),
    median_sample_zero_fraction = median(rowMeans(Xg == 0)),
    median_positive_count = median(rowSums(Xg > 0)),
    mean_positive_count = mean(rowSums(Xg > 0))
  )
}))
write_csv(group_summary, "group_sparsity_summary.csv")

X_ctr <- X[group == "CTR", , drop = FALSE]
X_crc <- X[group == "CRC", , drop = FALSE]
prev_ctr <- colMeans(X_ctr > 0)
prev_crc <- colMeans(X_crc > 0)
abs_diff <- abs(prev_crc - prev_ctr)

write_csv(
  data.frame(
    metric = c(
      "mean_abs_prevalence_difference",
      "median_abs_prevalence_difference",
      "max_abs_prevalence_difference",
      "features_abs_diff_gt_0.10",
      "features_abs_diff_gt_0.20",
      "features_abs_diff_gt_0.30"
    ),
    value = c(
      mean(abs_diff), median(abs_diff), max(abs_diff),
      sum(abs_diff > 0.10), sum(abs_diff > 0.20), sum(abs_diff > 0.30)
    )
  ),
  "prevalence_difference_summary.csv"
)

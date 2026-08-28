###############################################################################
## Prepare the SIAMCAT Zeller colorectal-cancer data used in the thesis
###############################################################################

source(file.path("applications", "zeller_crc", "config.R"))
source(file.path("applications", "zeller_crc", "helpers.R"))
check_zeller_root()

if (!requireNamespace("SIAMCAT", quietly = TRUE)) {
  stop(
    "Package 'SIAMCAT' is required. See the repository README for installation."
  )
}

data("feat_crc_zeller", package = "SIAMCAT")
data("meta_crc_zeller", package = "SIAMCAT")

get_loaded_object <- function(candidates) {
  for (name in candidates) {
    if (exists(name, inherits = TRUE)) return(get(name, inherits = TRUE))
  }
  stop("Could not find the Zeller data object in SIAMCAT.")
}

feat <- get_loaded_object(c("feat.crc.zeller", "feat_crc_zeller"))
meta <- get_loaded_object(c("meta.crc.zeller", "meta_crc_zeller"))
feat <- as.matrix(feat)
meta <- as.data.frame(meta)

if (!(group_column %in% colnames(meta))) {
  stop("Expected metadata column '", group_column, "'.")
}

# Align samples and retain only controls (CTR) and colorectal-cancer cases (CRC).
common_samples <- intersect(colnames(feat), rownames(meta))
if (length(common_samples) == 0L) stop("No matching sample IDs were found.")
feat <- feat[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]

keep_samples <- as.character(meta[[group_column]]) %in% group_order
feat <- feat[, keep_samples, drop = FALSE]
meta <- meta[keep_samples, , drop = FALSE]
group <- factor(as.character(meta[[group_column]]), levels = group_order)

sample_order <- order(group)
feat <- feat[, sample_order, drop = FALSE]
meta <- meta[sample_order, , drop = FALSE]
group <- droplevels(group[sample_order])

# scMANOVA expects samples x features.
X <- t(feat)
storage.mode(X) <- "numeric"
if (anyNA(X) || any(X < 0)) stop("Unexpected values in the Zeller feature matrix.")

original_features <- ncol(X)

# Remove the aggregate UNMAPPED feature.
if (remove_unmapped_feature) {
  is_unmapped <- grepl("^UNMAPPED$|^UNMAPPED\\b", colnames(X), ignore.case = TRUE)
} else {
  is_unmapped <- rep(FALSE, ncol(X))
}
removed_unmapped <- sum(is_unmapped)
X <- X[, !is_unmapped, drop = FALSE]
after_unmapped <- ncol(X)

# Remove features that are zero in every retained sample.
is_all_zero <- colSums(X > 0) == 0
removed_all_zero <- sum(is_all_zero)
X <- X[, !is_all_zero, drop = FALSE]
after_all_zero <- ncol(X)

# Thesis prevalence filter: retain features present in at least 5% of samples.
prevalence <- colMeans(X > 0)
keep_feature <- prevalence >= prevalence_cutoff
removed_by_prevalence <- sum(!keep_feature)
X <- X[, keep_feature, drop = FALSE]

n_vec <- as.numeric(table(group)[group_order])

prepared <- list(
  X = X,
  group = group,
  n = n_vec,
  group_order = group_order,
  metadata = meta,
  prevalence_cutoff = prevalence_cutoff
)

saveRDS(prepared, file.path(zeller_results_dir, "prepared_data.rds"))

write_csv(
  data.frame(group = names(table(group)), n = as.integer(table(group))),
  "group_counts.csv"
)

write_csv(
  data.frame(
    step = c(
      "original_features",
      "removed_unmapped",
      "after_unmapped_removal",
      "removed_all_zero",
      "after_all_zero_removal",
      "removed_by_prevalence_filter",
      "final_features"
    ),
    count = c(
      original_features,
      removed_unmapped,
      after_unmapped,
      removed_all_zero,
      after_all_zero,
      removed_by_prevalence,
      ncol(X)
    )
  ),
  "feature_filtering_counts.csv"
)

cat(
  "Prepared Zeller data:", nrow(X), "samples x", ncol(X), "features\n",
  "Group sizes:", paste(group_order, n_vec, sep = "=", collapse = ", "), "\n"
)

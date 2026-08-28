###############################################################################
## Zeller et al. colorectal-cancer application: thesis settings
###############################################################################

# Data preparation
prevalence_cutoff <- 0.05
remove_unmapped_feature <- TRUE

group_column <- "Group"
group_order <- c("CTR", "CRC")

# scMANOVA settings used for the thesis run
B_perm <- 199
lambda_bounds <- c(0, 100)
lambda0_bounds <- c(0, 100)
ident <- TRUE
p_value_perm <- TRUE
seed <- 20260710

# The analysis is intentionally single-core for portability.
ncpus <- 1
parallel_mode <- "no"

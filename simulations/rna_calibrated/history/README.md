# Archived RNA thesis results

`summary.csv` is the complete 1000-replication summary used for the RNA-calibrated thesis table. `scman_ind_h0.txt` and `scman_ind_h1.txt` are the combined raw scMAN-ind outputs retained as supporting historical results.

The file contents are unchanged; only the directory and file names were shortened to avoid Windows path-length problems.

From the repository root, recreate the exact thesis table with:

```r
source("simulations/rna_calibrated/make_thesis_table.R")
```

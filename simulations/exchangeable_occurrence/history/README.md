# Archived thesis results

This folder contains the retained historical outputs used for the three exchangeable beta-binomial occurrence simulation tables in the submitted thesis.

Only thesis cases are included: two groups, total `n = 10, 20`, `p = 50, 100, 150`, and the five H0/H1 scenarios. The method folders are `chen`, `scman_exc`, and `scman_ind`. Scenario folders/files use `h0`, `mean1`, `mean5`, `prob015`, and `prob030`.

Raw result contents and summary-table contents are unchanged; only folder and file names were shortened to avoid Windows path-length problems.

From the repository root, recreate the exact thesis tables with:

```r
source("simulations/exchangeable_occurrence/make_thesis_tables.R")
```

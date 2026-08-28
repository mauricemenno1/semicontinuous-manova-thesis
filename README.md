# Master's thesis replication repository

This repository contains the computational material used in the final master's thesis on tests for high-dimensional semicontinuous data. It was cleaned and audited against the submitted thesis source so that it contains the analyses, tables and computational figures that actually appear in the thesis, without obsolete exploratory branches.

It supports two complementary reproduction routes:

1. **Exact thesis tables from archived results.** The historical Monte Carlo outputs used for the submitted thesis are retained in each simulation study under `history/`. Use these when the goal is to reproduce the exact numerical values printed in the thesis.
2. **Fresh end-to-end reruns.** Clean, single-core scripts rerun the simulations from scratch. They replace the original split/multicore workflows and are intended to be portable to Windows.

The repository covers all **7 computational tables** and **10 computational figures** in the final thesis.

## Repository structure

```text
.
├── R/                              Shared fitting, simulation and table helpers
├── applications/
│   └── zeller_crc/                 Zeller colorectal-cancer application
├── data/
│   └── rna_simulation_parameters.RData
├── figures/                        Standalone thesis figure scripts
├── package/                        Edited Bernoulli-enabled semicontMANOVA package
├── simulations/
│   ├── independent_occurrence/     Independent-Bernoulli injection study
│   ├── exchangeable_occurrence/    Exchangeable beta-binomial injection study
│   └── rna_calibrated/             RNA-calibrated simulation study
├── make_thesis_tables.R            Rebuild all seven tables from archived results
└── README.md
```

The historical result hierarchy and raw filenames were deliberately shortened to avoid the traditional Windows path-length limit. The numerical contents of the retained historical result files were not changed by this path cleanup.

## Exact thesis tables from archived results

From the repository root, run:

```r
source("make_thesis_tables.R")
```

No Monte Carlo simulations are rerun. The script reads the archived thesis outputs and writes all seven tables as both `.csv` and `.tex` files.

Outputs are written to:

```text
simulations/independent_occurrence/tables/thesis/
├── size.csv
├── size.tex
├── mean_power.csv
├── mean_power.tex
├── occurrence_power.csv
└── occurrence_power.tex

simulations/exchangeable_occurrence/tables/thesis/
├── size.csv
├── size.tex
├── mean_power.csv
├── mean_power.tex
├── occurrence_power.csv
└── occurrence_power.tex

simulations/rna_calibrated/tables/thesis/
├── rna_simulation.csv
└── rna_simulation.tex
```

These values were checked against the submitted thesis tables.

### Historical-result layout

For the two main simulation studies, archived results are stored as:

```text
history/
├── chen/
├── scman_exc/
└── scman_ind/
```

Each method folder contains five summary CSVs and corresponding raw-result folders:

- `h0`: null hypothesis;
- `mean1`: mean shift with `c1 = 1`;
- `mean5`: mean shift with `c1 = 5`;
- `prob015`: occurrence-probability shift with `c2 = 0.15`;
- `prob030`: occurrence-probability shift with `c2 = 0.30`.

Only final-thesis cases are retained: two groups, total sample size `n = 10, 20` (`n_k = 5, 10` per group), `p = 50, 100, 150`, and the valid parameter combinations displayed in the thesis. Obsolete `p = 200`, `K = 4`, recovery/progress material and old split-run scripts were removed.

One archived independent-occurrence summary value is `NA` even though its corresponding raw 250-replication `.txt` file is complete. The table reader uses that raw file only for this one missing summary entry and recovers the thesis value `1.000`.

For the RNA study, `simulations/rna_calibrated/history/summary.csv` is the complete 1000-replication summary used for the thesis table. The two combined raw scMAN-ind files are retained as supporting historical output.

## Installation

Run R or RStudio with the repository root as the working directory.

Install the edited thesis version of `semicontMANOVA` with:

```r
source("package/install_package.R")
```

The script installs the included source package and checks that `scMANOVA()` contains the Bernoulli-enabled `missing.model` argument. It also installs the required CRAN dependencies `matrixcalc` and `mvtnorm` when missing.

On Windows, installing an R package from source may require a compatible version of Rtools.

The Zeller application additionally requires the Bioconductor package `SIAMCAT`, because the analysis loads the public `feat_crc_zeller` and `meta_crc_zeller` datasets from that package.

## Fresh simulation reruns

The fresh runners are intentionally **single-core**. They process the full thesis grid serially and write restartable output under each study's `results/raw/` directory. Newly generated raw results are ignored by Git by default.

### 1. Independent-Bernoulli occurrence simulation

```r
source("simulations/independent_occurrence/run_all.R")
```

Thesis settings:

- two groups;
- `n_k = 5, 10` per group;
- `p = 50, 100, 150`;
- baseline occurrence/zero-probability settings `0.2, 0.5, 0.8`;
- Gaussian correlation `rho = 0, 0.4`;
- mean shifts `c1 = 1, 5`;
- occurrence-probability shifts `c2 = 0.15, 0.30`, omitting invalid probability combinations;
- `250` Monte Carlo replications;
- `B = 100` permutations.

The same generated dataset within each replication is supplied to Chen, scMAN-exc and scMAN-ind. The runner processes all `p = 50` cases first, then `p = 100`, then `p = 150`, and produces tables for empirical size, mean-shift power and occurrence-shift power.

### 2. Exchangeable beta-binomial occurrence simulation

```r
source("simulations/exchangeable_occurrence/run_all.R")
```

This uses the same thesis grid and the same three methods, but simulates the occurrence component using the exchangeable beta-binomial mechanism from the thesis. The occurrence ICC is fixed at `0.05`.

It produces the same three table types as the independent-occurrence study.

### 3. RNA-calibrated simulation

```r
source("simulations/rna_calibrated/run_all.R")
```

Thesis settings:

- two groups;
- parameters loaded from `data/rna_simulation_parameters.RData`;
- `1000` Monte Carlo replications;
- `B = 1000` permutations;
- lambda search intervals `[0, 100]` under both `H1` and `H0`.

The output is the single RNA H0/H1 comparison table used in the thesis.

### Reproducibility note

The clean runners use fixed seeds and reproduce the same simulation design. The original thesis computations were sometimes executed in separate parts and merged afterward, so a fresh serial rerun can consume random numbers in a different order or differ slightly because of R/package/platform numerical behavior.

Therefore:

- use the archived `history/` results for the **exact submitted thesis values**;
- use the clean runners to **reproduce the analyses from scratch**.

## Zeller colorectal-cancer application

Run:

```r
source("applications/zeller_crc/run_all.R")
```

The workflow:

1. loads the Zeller colorectal-cancer data from `SIAMCAT`;
2. retains CTR and CRC samples;
3. removes `UNMAPPED` and all-zero features;
4. applies the 5% prevalence filter used in the thesis;
5. fits scMAN-exc and scMAN-ind with the thesis settings;
6. writes descriptive/model summaries;
7. regenerates the four Zeller figures used in the thesis.

The prepared thesis dataset contains 141 samples (88 CTR and 53 CRC) and 368 retained features. The Zeller application does not produce a thesis table; its results appear in the text and figures.

Generated figures are written to `applications/zeller_crc/figures/` and numerical summaries to `applications/zeller_crc/results/`.

## Standalone thesis figures

The remaining computational figures are generated with:

```r
source("figures/01_count_distribution_comparison.R")
source("figures/02_occurrence_fit_illustration.R")
source("figures/03_power_loss_margin.R")
```

Together with the four Zeller figures, these scripts account for all 10 computational figures in the final thesis. Reference copies of the final thesis figures are retained in `figures/generated/` and `applications/zeller_crc/figures/`.

## Thesis-output map

| Thesis output | Reproduction source |
| --- | --- |
| Independent occurrence: empirical size | `simulations/independent_occurrence/make_thesis_tables.R` |
| Independent occurrence: mean-shift power | same |
| Independent occurrence: occurrence-shift power | same |
| Exchangeable occurrence: empirical size | `simulations/exchangeable_occurrence/make_thesis_tables.R` |
| Exchangeable occurrence: mean-shift power | same |
| Exchangeable occurrence: occurrence-shift power | same |
| RNA-calibrated H0/H1 table | `simulations/rna_calibrated/make_thesis_table.R` |
| Binomial vs beta-binomial count distribution | `figures/01_count_distribution_comparison.R` |
| Four occurrence-model fit illustrations | `figures/02_occurrence_fit_illustration.R` |
| Power-loss margin figure | `figures/03_power_loss_margin.R` |
| Four Zeller CRC figures | `applications/zeller_crc/make_figures.R` via `run_all.R` |

## Suggested reproduction workflow

For a reader who only wants the submitted thesis outputs:

```r
# Exact numerical tables from archived thesis runs
source("make_thesis_tables.R")

# Zeller application and its four figures
source("applications/zeller_crc/run_all.R")

# Remaining standalone thesis figures
source("figures/01_count_distribution_comparison.R")
source("figures/02_occurrence_fit_illustration.R")
source("figures/03_power_loss_margin.R")
```

For a full fresh rerun, install the package first and then run the three simulation `run_all.R` scripts described above.

## Windows and ZIP extraction

This repository intentionally uses compact archived-result paths. The longest path inside the release ZIP is kept well below the traditional Windows `MAX_PATH` limit, leaving room for extraction under a normal directory such as `C:\Users\<name>\Downloads\`.

If you move the repository inside an unusually deep directory tree, shorter locations such as `C:\thesis-replication\` are still preferable for R projects in general.

## Generated files and Git

`.gitignore` excludes restartable outputs from **new** Monte Carlo runs, local R/RStudio files and reproducible intermediate data. The thesis `history/` archives and retained final tables/figures are intended to be committed to GitHub.

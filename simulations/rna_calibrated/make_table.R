###############################################################################
## Make the RNA table used in the thesis
###############################################################################
source(file.path("simulations", "rna_calibrated", "config.R"))
assert_repository_root()
if (!file.exists(SUMMARY_FILE)) stop("Run summarize_results.R first.")

tab <- read.csv(SUMMARY_FILE, stringsAsFactors = FALSE)
model_order <- c("scMAN_MV", "scMAN_Bernoulli")
scenario_order <- c("H0", "H1")
tab$model_rank <- match(tab$model, model_order)
tab$scenario_rank <- match(tab$scenario, scenario_order)
tab <- tab[order(tab$model_rank, tab$scenario_rank), , drop = FALSE]

model_csv <- c(scMAN_MV = "scMAN-exc", scMAN_Bernoulli = "scMAN-ind")
model_tex <- c(
  scMAN_MV = "$\\mathrm{scMAN}_{\\mathrm{exc}}$",
  scMAN_Bernoulli = "$\\mathrm{scMAN}_{\\mathrm{ind}}$"
)

out_csv <- data.frame(
  Method = unname(model_csv[tab$model]),
  Scenario = tab$scenario,
  n = tab$n,
  p = tab$p,
  rejection_proportion = tab$rejection_proportion,
  mean_lambda = tab$mean_lambda,
  mean_lambda0 = tab$mean_lambda0,
  stringsAsFactors = FALSE
)

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(out_csv, file.path(TABLE_DIR, "rna_simulation.csv"), row.names = FALSE, na = "NA")

fmt3 <- function(x) ifelse(is.na(x), "x", sprintf("%.3f", x))
fmt2 <- function(x) ifelse(is.na(x), "x", sprintf("%.2f", x))

lines <- c(
  "\\begin{table}[H]",
  "\\scriptsize",
  "\\centering",
  "\\caption{RNA simulation results for $\\mathrm{scMAN}_{\\mathrm{exc}}$ and $\\mathrm{scMAN}_{\\mathrm{ind}}$.}",
  "\\label{tab:rna_dataset}",
  "\\begin{tabular}{llrrrrr}",
  "\\toprule",
  "Method & Scenario & $n$ & $p$ & Rej. prop. & $\\bar\\lambda$ & $\\bar\\lambda_0$ \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(tab))) {
  lines <- c(lines, paste0(
    unname(model_tex[tab$model[i]]), " & $H_", ifelse(tab$scenario[i] == "H0", "0", "1"), "$ & ",
    tab$n[i], " & ", tab$p[i], " & ", fmt3(tab$rejection_proportion[i]), " & ",
    fmt2(tab$mean_lambda[i]), " & ", fmt2(tab$mean_lambda0[i]), " \\\\" 
  ))
}

lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(lines, file.path(TABLE_DIR, "rna_simulation.tex"))
cat("Wrote thesis RNA table to:", TABLE_DIR, "\n")

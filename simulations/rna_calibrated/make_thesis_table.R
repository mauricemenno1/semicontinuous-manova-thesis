###############################################################################
## Recreate the submitted-thesis RNA table from the archived historical output
###############################################################################

historical_root <- file.path("simulations", "rna_calibrated", "history")
summary_path <- file.path(historical_root, "summary.csv")
out_dir <- file.path("simulations", "rna_calibrated", "tables", "thesis")

if (!file.exists(summary_path)) stop("Historical RNA summary not found: ", summary_path)

x <- read.csv(summary_path, stringsAsFactors = FALSE)

## The archived summary contains all 1000 completed replications for both H0 and H1 for both methods.
if (!all(x$n_simEnded == 1000L)) {
  stop("Historical RNA summary is not complete: expected 1000 replications per row.")
}

model_key <- ifelse(
  x$Model == "Original multiv. Bernoulli", "scMAN_MV",
  ifelse(grepl("Independent Bernoulli", x$Model, fixed = TRUE), "scMAN_Bernoulli", NA_character_)
)
if (anyNA(model_key)) stop("Unexpected model name in historical RNA summary.")

scenario <- ifelse(x$Scenario == "$H_0$", "H0", ifelse(x$Scenario == "$H_1$", "H1", NA_character_))
if (anyNA(scenario)) stop("Unexpected scenario name in historical RNA summary.")

tab <- data.frame(
  model = model_key,
  scenario = scenario,
  n = as.integer(x$n),
  p = as.integer(x$p),
  rejection_proportion = as.numeric(x$rejPERM),
  mean_lambda = as.numeric(x$meanLambda),
  mean_lambda0 = as.numeric(x$meanLambda0),
  n_completed = as.integer(x$n_simEnded),
  stringsAsFactors = FALSE
)

model_order <- c("scMAN_MV", "scMAN_Bernoulli")
scenario_order <- c("H0", "H1")
tab <- tab[order(match(tab$model, model_order), match(tab$scenario, scenario_order)), , drop = FALSE]

model_csv <- c(scMAN_MV = "scMAN-exc", scMAN_Bernoulli = "scMAN-ind")
model_tex <- c(
  scMAN_MV = "$\\mathrm{scMAN}_{\\mathrm{exc}}$",
  scMAN_Bernoulli = "$\\mathrm{scMAN}_{\\mathrm{ind}}$"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
write.csv(out_csv, file.path(out_dir, "rna_simulation.csv"), row.names = FALSE, na = "NA")

fmt3 <- function(v) sprintf("%.3f", v)
fmt2 <- function(v) sprintf("%.2f", v)

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
writeLines(lines, file.path(out_dir, "rna_simulation.tex"))

cat("Wrote historical RNA thesis table to:", out_dir, "\n")

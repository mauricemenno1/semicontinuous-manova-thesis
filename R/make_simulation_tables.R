###############################################################################
## Build the three simulation tables used in the thesis
###############################################################################

simulation_table_scenarios <- function(table_id) {
  if (table_id == 1L) {
    out <- expand.grid(
      c1 = 0,
      c2 = 0,
      pmiss = c(0.2, 0.5, 0.8),
      rho = c(0, 0.4),
      KEEP.OUT.ATTRS = FALSE
    )
  } else if (table_id == 2L) {
    out <- expand.grid(
      c1 = c(1, 5),
      c2 = 0,
      pmiss = c(0.2, 0.5, 0.8),
      rho = c(0, 0.4),
      KEEP.OUT.ATTRS = FALSE
    )
  } else if (table_id == 3L) {
    out <- rbind(
      expand.grid(c1 = 0, c2 = 0.15, pmiss = c(0.2, 0.5, 0.8), rho = c(0, 0.4), KEEP.OUT.ATTRS = FALSE),
      expand.grid(c1 = 0, c2 = 0.30, pmiss = c(0.2, 0.5), rho = c(0, 0.4), KEEP.OUT.ATTRS = FALSE)
    )
  } else {
    stop("table_id must be 1, 2, or 3.")
  }
  out[order(out$c1, out$c2, out$pmiss, out$rho), , drop = FALSE]
}

## Thesis order: exchangeable scMANOVA, independent scMANOVA, Chen.
model_order <- c("scMAN_MV", "scMAN_Bernoulli", "CHEN")
model_csv_names <- c(
  scMAN_MV = "scMAN-exc",
  scMAN_Bernoulli = "scMAN-ind",
  CHEN = "Chen"
)
model_latex_names <- c(
  scMAN_MV = "$\\mathrm{scMAN}_{\\mathrm{exc}}$",
  scMAN_Bernoulli = "$\\mathrm{scMAN}_{\\mathrm{ind}}$",
  CHEN = "Chen"
)

lookup_rejection <- function(summary_df, c1, c2, pmiss, rho, model, nk, p) {
  z <- summary_df[
    abs(summary_df$c1 - c1) < 1e-12 &
      abs(summary_df$c2 - c2) < 1e-12 &
      abs(summary_df$pmiss - pmiss) < 1e-12 &
      abs(summary_df$rho - rho) < 1e-12 &
      summary_df$model == model &
      summary_df$n_per_group == nk &
      summary_df$p == p,
    , drop = FALSE
  ]
  if (nrow(z) == 0L) return(NA_real_)
  z$rejection_proportion[1]
}

build_simulation_table <- function(summary_df, table_id) {
  settings <- simulation_table_scenarios(table_id)
  available_models <- model_order[model_order %in% unique(summary_df$model)]
  nks <- c(5L, 10L)
  ps <- c(50L, 100L, 150L)
  rows <- list()

  for (i in seq_len(nrow(settings))) {
    setting <- settings[i, , drop = FALSE]
    for (model in available_models) {
      row <- data.frame(
        c1 = setting$c1,
        c2 = setting$c2,
        pmiss = setting$pmiss,
        rho = setting$rho,
        model_key = model,
        model = unname(model_csv_names[model]),
        stringsAsFactors = FALSE
      )
      for (nk in nks) {
        for (p in ps) {
          col <- paste0("nk", nk, "_p", p)
          row[[col]] <- lookup_rejection(
            summary_df, setting$c1, setting$c2, setting$pmiss,
            setting$rho, model, nk, p
          )
        }
      }
      rows[[length(rows) + 1L]] <- row
    }
  }

  do.call(rbind, rows)
}

fmt_value <- function(x) ifelse(is.na(x), "x", sprintf("%.3f", x))
fmt_c1 <- function(x) sprintf("%.0f", x)
fmt_c2 <- function(x) {
  if (abs(x) < 1e-12) return("0")
  if (abs(x - 0.30) < 1e-12) return("0.30")
  sprintf("%.2f", x)
}
fmt_one_decimal_or_zero <- function(x) if (abs(x) < 1e-12) "0" else sprintf("%.1f", x)

simulation_table_metadata <- function(study_key, table_id) {
  if (!study_key %in% c("independent", "exchangeable")) {
    stop("study_key must be 'independent' or 'exchangeable'.")
  }

  if (study_key == "independent") {
    captions <- c(
      "Empirical rejection proportions under the null hypothesis with independent Bernoulli occurrence generation. The nominal significance level is $0.05$.",
      "Empirical power under differences in the continuous means with independent Bernoulli occurrence generation.",
      "Empirical power under differences in occurrence probabilities with independent Bernoulli occurrence generation."
    )
    labels <- c(
      "tab:size-independent",
      "tab:mean-power-independent",
      "tab:occurrence-power-independent"
    )
  } else {
    captions <- c(
      "Empirical rejection proportions under the null hypothesis with exchangeable beta-binomial occurrence generation. The nominal significance level is $0.05$.",
      "Empirical power under differences in the continuous means with exchangeable beta-binomial occurrence generation.",
      "Empirical power under differences in occurrence probabilities with exchangeable beta-binomial occurrence generation. Results are calculated over the completed replications in each scenario."
    )
    labels <- c(
      "tab:size-exchangeable",
      "tab:mean-power-exchangeable",
      "tab:occurrence-power-exchangeable"
    )
  }

  filenames <- c("size", "mean_power", "occurrence_power")
  list(caption = captions[table_id], label = labels[table_id], filename = filenames[table_id])
}

write_table_latex <- function(tab, path, caption, label) {
  value_cols <- c(
    "nk5_p50", "nk5_p100", "nk5_p150",
    "nk10_p50", "nk10_p100", "nk10_p150"
  )

  lines <- c(
    "\\begin{table}[H]",
    "\\scriptsize",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular}{lllllllllll}",
    "\\hline",
    "$c_1$ & $c_2$ & $\\pi_{j1}$ & $\\gamma$ & Method & $p=50$ & $p=100$ & $p=150$ & $p=50$ & $p=100$ & $p=150$ \\\\",
    "\\hline",
    "& & & & & $n_k=5$ & $n_k=5$ & $n_k=5$ & $n_k=10$ & $n_k=10$ & $n_k=10$ \\\\" 
  )

  previous_setting <- NULL
  for (i in seq_len(nrow(tab))) {
    current_setting <- paste(tab$c1[i], tab$c2[i], tab$pmiss[i], tab$rho[i], sep = "|")
    show_setting <- is.null(previous_setting) || current_setting != previous_setting

    setting_cells <- if (show_setting) {
      c(
        fmt_c1(tab$c1[i]),
        fmt_c2(tab$c2[i]),
        fmt_one_decimal_or_zero(tab$pmiss[i]),
        fmt_one_decimal_or_zero(tab$rho[i])
      )
    } else {
      rep("", 4)
    }

    vals <- vapply(value_cols, function(v) fmt_value(tab[[v]][i]), character(1))
    method <- unname(model_latex_names[tab$model_key[i]])

    lines <- c(
      lines,
      paste0(
        paste(setting_cells, collapse = " & "), " & ", method, " & ",
        paste(vals, collapse = " & "), " \\\\" 
      )
    )
    previous_setting <- current_setting
  }

  lines <- c(lines, "\\hline", "\\end{tabular}", "\\end{table}")
  writeLines(lines, path)
}

write_simulation_tables <- function(summary_df, table_dir, study_key) {
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  for (table_id in 1:3) {
    meta <- simulation_table_metadata(study_key, table_id)
    tab <- build_simulation_table(summary_df, table_id)

    csv_tab <- tab[, setdiff(names(tab), "model_key"), drop = FALSE]
    write.csv(
      csv_tab,
      file.path(table_dir, paste0(meta$filename, ".csv")),
      row.names = FALSE,
      na = "NA"
    )

    write_table_latex(
      tab,
      file.path(table_dir, paste0(meta$filename, ".tex")),
      meta$caption,
      meta$label
    )
  }

  invisible(TRUE)
}

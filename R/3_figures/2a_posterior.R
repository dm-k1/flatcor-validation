## =====================================================================
##  Figures of study 2a, from its results file. Nothing is computed here.
##    Rscript R/3_figures/2a_posterior.R
## =====================================================================
library(ggplot2)
source("R/shared/figures.R") # read_results(), figure_folder

results <- read_results("2a")
stats <- as.data.frame(results$tables$sbc_stats)


## The fits numbered through per prior and setting: the SBC package numbers
## them (sim_id) within each batch, so in a run of several batches the same
## sim_id occurs once per batch
stats$sim_id <- ave(
  seq_len(nrow(stats)),
  stats$prior, stats$k, stats$variable,
  FUN = seq_along
)


## Simulation-based calibration: the ECDF difference of the ranks (the SBC
## package's plot), one panel per prior and setting: a column per setting,
## a row per prior
panels <- list()
for (setting_k in sort(unique(stats$k))) {
  for (prior_name in unique(stats$prior)) {
    cell <- stats[stats$prior == prior_name & stats$k == setting_k, ]
    panels[[length(panels) + 1]] <- SBC::plot_ecdf_diff(cell) +
      labs(title = sprintf("%s, k = %d, n = %d", prior_name, setting_k, cell$n[1]))
  }
}

figure <- patchwork::wrap_plots(panels, ncol = length(unique(stats$k)), byrow = FALSE) +
  patchwork::plot_annotation(
    title = "Simulation-based calibration: ECDF difference of the ranks"
  )

ggsave(
  file.path(figure_folder, "2a_sbc_ecdf_difference.png"), figure,
  width = 7 * length(unique(stats$k)), height = 14, dpi = 150
)

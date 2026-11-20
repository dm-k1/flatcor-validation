## =====================================================================
##  Figures of study 2a, from its per-dataset file. Nothing is computed here.
##    Rscript R/3_figures/2a_posterior.R
## =====================================================================
library(ggplot2)
source("R/shared/figures.R") # read_per_dataset_file(), figure_folder

stats <- read_per_dataset_file("2a_sbc_fits") # one row per variable of each SBC fit


## Each fit by its dataset: the SBC package numbers fits (sim_id) within
## each batch, the replication numbers them through
stats$sim_id <- stats$replication


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

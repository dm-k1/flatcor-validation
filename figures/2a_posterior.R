## =====================================================================
##  Figures of study 2a, from its results file. Nothing is computed here.
##    Rscript figures/2a_posterior.R
## =====================================================================
library(ggplot2)
source("R/shared/figures.R") # read_results(), figure_folder

results <- read_results("2a")
stats <- results$tables$sbc_stats


## Simulation-based calibration: the ECDF difference of the ranks (the SBC
## package's plot), one panel row per prior, for the first setting
first_k <- stats$k[1]

panels <- lapply(unique(stats$prior), function(prior_name) {
  of_prior <- stats[stats$prior == prior_name & stats$k == first_k, ]
  SBC::plot_ecdf_diff(as.data.frame(of_prior)) + labs(title = prior_name)
})

figure <- patchwork::wrap_plots(panels, ncol = 1) +
  patchwork::plot_annotation(
    title = sprintf("Simulation-based calibration, k = %d: ECDF difference of the ranks", first_k)
  )

ggsave(
  file.path(figure_folder, "2a_sbc_ecdf_difference.png"), figure,
  width = 10, height = 14, dpi = 150
)

## =====================================================================
##  Figures of study 3a, from its results file. Nothing is computed here.
##    Rscript R/3_figures/3a_bayes_factors.R
##  What the figure shows and why: docs/validation_tests.md, study 3a.
## =====================================================================
library(ggplot2)
source("R/shared/figures.R") # read_results(), figure_folder

results <- read_results("3a")
bins <- results$tables$calibration_bins
if (nrow(bins) == 0) stop("No bin has enough datasets: nothing to draw.")


## Calibration: per bin of log BF10, the observed share of datasets without
## x1 against the expected share, 1 / (1 + BF10). On the diagonal if the
## Bayes factors are right. Grey: 2 Monte Carlo SD under calibration.
figure <- ggplot(bins, aes(expected, estimate)) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_errorbar(
    aes(ymin = expected - 2 * mcse, ymax = expected + 2 * mcse),
    width = 0,
    colour = "grey70"
  ) +
  geom_point(size = 2, colour = "#0072B2") +
  geom_line(colour = "#0072B2") +
  facet_grid(prior ~ cell) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "expected share without x1: mean of 1 / (1 + BF10) in the bin",
    y = "observed share without x1",
    title = "Calibration of the Bayes factor for x1",
    subtitle = "on the diagonal if the Bayes factors are right; grey: 2 Monte Carlo SD"
  ) +
  theme_bw()

ggsave(
  file.path(figure_folder, "3a_bayes_factor_calibration.png"), figure,
  width = 9, height = 9, dpi = 150
)

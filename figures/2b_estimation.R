## =====================================================================
##  Figures of study 2b, from its results file. Nothing is computed here.
##    Rscript figures/2b_estimation.R
## =====================================================================
library(ggplot2)
source("R/shared/figures.R") # read_results(), figure_folder
source("R/shared/comparison_methods.R") # method_scales(): the colours of the methods

results <- read_results("2b")
table <- results$tables$estimates


## Each estimate against its reference (dotted), per cell, with 2 Monte Carlo SE
figure_of <- function(quantity, ylab, title, file) {
  rows <- table[table$quantity == quantity, ]

  plot <- ggplot(rows, aes(n, estimate, colour = method)) +
    geom_hline(aes(yintercept = expected), linetype = "dotted") +
    geom_pointrange(
      aes(ymin = estimate - 2 * mcse, ymax = estimate + 2 * mcse),
      position = position_dodge(width = 0.15),
      size = 0.25,
      alpha = 0.8
    ) +
    scale_x_log10(breaks = c(50, 200)) +
    method_scales(rows$method) +
    facet_grid(
      others + rho ~ k,
      labeller = labeller(
        k = function(v) paste("k =", v),
        rho = function(v) paste("correlation =", v)
      )
    ) +
    labs(
      x = "sample size (log scale)",
      y = ylab,
      title = title,
      subtitle = "bars: 2 Monte Carlo SE; dotted: the reference"
    ) +
    theme_bw()

  ggsave(file.path(figure_folder, file), plot, width = 9, height = 10, dpi = 150)
}


figure_of(
  "coverage of x1's 95% interval",
  ylab = "coverage of x1's 95% interval",
  title = "Do the intervals contain the true slope of x1 (0.2)?",
  file = "2b_coverage.png"
)
figure_of(
  "bias of x1 (share of 0.2)",
  ylab = "bias of x1, as a share of its true value",
  title = "Bias: how far estimates of x1 are pulled towards zero",
  file = "2b_bias.png"
)
figure_of(
  "prediction error / OLS",
  ylab = "excess prediction error / OLS",
  title = "Do the posterior means predict new data better than OLS?",
  file = "2b_prediction.png"
)

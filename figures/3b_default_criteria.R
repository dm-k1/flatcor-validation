## =====================================================================
##  Figures of study 3b, from its results file. Nothing is computed here.
##    Rscript figures/3b_default_criteria.R
## =====================================================================
library(ggplot2)
source("R/shared/figures.R") # read_results(), figure_folder
source("R/shared/comparison_methods.R") # method_scales(): the colours of the methods

results <- read_results("3b")
medians <- results$tables$growing_medians


## Growing k (3b.8, 3b.9): the median log Bayes factor for the truth against
## n; expected to grow for every method
figure <- ggplot(medians, aes(n, median_log_bf, colour = method)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  geom_point(size = 1.3, alpha = 0.8) +
  scale_x_log10(breaks = sort(unique(medians$n))) +
  method_scales(medians$method) +
  facet_grid(
    truth ~ growth,
    labeller = labeller(truth = c(null = "null true", true = "3-predictor model true"))
  ) +
  labs(
    x = "n (log scale)",
    y = "median log Bayes factor for the truth",
    title = "Consistency when the number of predictors grows with n (k at most 90)"
  ) +
  theme_bw()

ggsave(file.path(figure_folder, "3b_growing_k.png"), figure, width = 11, height = 7, dpi = 150)

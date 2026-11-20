## =====================================================================
##  Figures of study 1b, from its results file. Nothing is computed here.
##    Rscript R/3_figures/1b_invariance.R
## =====================================================================
library(ggplot2)
source("R/shared/figures.R") # read_results(), figure_folder
source("R/shared/comparison_methods.R") # contenders, method_scales(): names and colours

results <- read_results("1b")
pooled <- as.data.frame(results$tables$pooled)
pooled$prior <- names(contenders)[match(pooled$prior, contenders)]


## The mean change in log Bayes factor per prior and change of units or
## coding, with 2 Monte Carlo SEs either side, inside the band of changes
## that would matter (plus or minus delta, docs/validation_tests.md)
figure <- ggplot(pooled, aes(estimate, quantity, colour = prior)) +
  annotate(
    "rect", xmin = -pooled$delta[1], xmax = pooled$delta[1], ymin = -Inf, ymax = Inf,
    fill = "grey90"
  ) +
  geom_vline(xintercept = 0) +
  geom_pointrange(
    aes(xmin = estimate - 2 * mcse, xmax = estimate + 2 * mcse),
    position = position_dodge(width = 0.6)
  ) +
  method_scales(pooled$prior) +
  facet_wrap(~bayes_factor, ncol = 1, scales = "free_y") +
  labs(
    x = paste(
      "mean change in log Bayes factor",
      "(2 Monte Carlo SEs either side; grey: a change that matters)"
    ),
    y = NULL,
    title = "Invariance: do units and coding change the Bayes factor?"
  ) +
  theme_bw()

ggsave(file.path(figure_folder, "1b_pooled_changes.png"), figure, width = 10, height = 6, dpi = 150)

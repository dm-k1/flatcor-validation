## =====================================================================
##  Figures of study 1a, from its results file. Nothing is computed here.
##    Rscript R/3_figures/1a_prior_sampling.R
## =====================================================================
library(ggplot2)
source("R/shared/figures.R") # read_results(), figure_folder
source("R/shared/comparison_methods.R") # contenders, method_scales(): names and colours

results <- read_results("1a")
quantiles <- as.data.frame(results$tables$quantiles)
quantiles$prior <- names(contenders)[match(quantiles$prior, contenders)]


## Every percentile of the prior draws: its distance from the expected value
## in Monte Carlo SEs, one panel per quantity; if the model draws the prior,
## about 1 in 20 lies beyond 2 and almost none beyond 3
figure <- ggplot(quantiles, aes(factor(percentile), distance_in_mcse, colour = prior)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = c(-3, 3), linetype = "dotted", colour = "grey50") +
  geom_point(position = position_dodge(width = 0.6), size = 1.8) +
  method_scales(quantiles$prior) +
  facet_wrap(~ paste0(quantity, " (", design, ")"), ncol = 2) +
  labs(
    x = "percentile",
    y = "estimate minus expected value, in Monte Carlo SEs",
    title = "Prior draws against the exact prior: every percentile within its Monte Carlo error?"
  ) +
  theme_bw()

ggsave(
  file.path(figure_folder, "1a_percentiles_in_mcse.png"), figure,
  width = 10, height = 9, dpi = 150
)

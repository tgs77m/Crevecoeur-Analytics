# ---------------------------------------------------------
# Libraries
# ---------------------------------------------------------
library(ggplot2)
library(ggridges)
library(cowplot)
library(dplyr)
library(tidyr)

# ---------------------------------------------------------
# Exact function
# ---------------------------------------------------------
exacts <- function(x) (x^2 + 2*x)^0.9

# ---------------------------------------------------------
# Generalized Binomial Expansions (large-x branch)
# ---------------------------------------------------------

# 2nd-order GBE
order2 <- function(x) {
  x^1.8 +
    1.8000 * x^0.8 -
    0.1800 * x^-0.2
}

# 3rd-order GBE
order3 <- function(x) {
  x^1.8 +
    1.8000 * x^0.8 -
    0.1800 * x^-0.2 +
    0.0120 * x^-1.2
}

# 4th-order GBE
order4 <- function(x) {
  x^1.8 +
    1.8000 * x^0.8 -
    0.1800 * x^-0.2 +
    0.0120 * x^-1.2 -
    0.00072 * x^-2.2
}

# ---------------------------------------------------------
# Corrected GBE2 (regression-based global correction)
# ---------------------------------------------------------
gbe2_corr <- function(x) {
  x^1.8 +
    1.8000 * x^0.8 -
    0.1800 * x^-0.2 -
    (4.907311e-06 * x - 1.175748e-03)
}

# ---------------------------------------------------------
# Domain and sampling
# ---------------------------------------------------------
set.seed(789)
x <- seq(40, 240, 1)
x100 <- sample(x, size = 100, replace = FALSE)

# ---------------------------------------------------------
# Compute approximations
# ---------------------------------------------------------
gbe2_raw <- order2(x100)
gbe3_raw <- order3(x100)
gbe4_raw <- order4(x100)
gbe2_corr_raw <- gbe2_corr(x100)

# ---------------------------------------------------------
# Signed errors (biases)
# ---------------------------------------------------------
gbe2_bias <- gbe2_raw - exacts(x100)
gbe3_bias <- gbe3_raw - exacts(x100)
gbe4_bias <- gbe4_raw - exacts(x100)
gbe2_corr_bias <- gbe2_corr_raw - exacts(x100)

# ---------------------------------------------------------
# Assemble long-format data for plotting + Wilcoxon
# ---------------------------------------------------------
bias_compare <- data.frame(
  x = x100,
  gbe2_corr = gbe2_corr_bias,
  gbe3 = gbe3_bias,
  gbe4 = gbe4_bias
) %>%
  pivot_longer(cols = c("gbe2_corr", "gbe3", "gbe4"),
               names_to = "condition",
               values_to = "value")

# ---------------------------------------------------------
# Ridge plot
# ---------------------------------------------------------
ggplot(bias_compare, aes(x = value, y = condition, fill = condition)) +
  geom_density_ridges(alpha = 0.6) +
  theme_cowplot() +
  ggtitle("Signed Error Comparison: Corrected GBE2 vs GBE3 vs GBE4") +
  xlab("Signed Error") +
  ylab("Ridge height = Local Density") +
  geom_vline(xintercept = 0, linetype = "dashed")

# ---------------------------------------------------------
# Boxplot comparison
# ---------------------------------------------------------
ggplot(bias_compare, aes(x = condition, y = value, fill = condition)) +
  geom_boxplot(notch = TRUE) +
  theme_cowplot() +
  ggtitle("Bias Comparison: Corrected GBE2 vs GBE3 vs GBE4") +
  geom_hline(yintercept = 0, linetype = "dashed")

# ---------------------------------------------------------
# Pairwise Wilcoxon tests (Holm correction)
# ---------------------------------------------------------
pairwise.wilcox.test(
  x = bias_compare$value,
  g = bias_compare$condition,
  paired = TRUE,
  p.adjust.method = "holm"
)

# ---------------------------------------------------------
# Friedman test for repeated-measures comparison
# ---------------------------------------------------------

# reshape to wide format for friedman.test
friedman_df <- bias_compare %>%
  select(x, condition, value) %>%
  pivot_wider(names_from = condition, values_from = value) %>%
  arrange(x)

# run Friedman test
friedman.test(
  y = as.matrix(friedman_df[, c("gbe2_corr", "gbe3", "gbe4")])
)


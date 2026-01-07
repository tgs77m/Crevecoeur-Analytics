# rm(list = ls())

library(ggridges)
library(ggplot2)
library(cowplot)

# 2nd‑order generalized binomial expansion (large‑x branch, valid for x > 2)
order2 <- function(x) x^1.8 +
  1.8000 * x^0.8 -
  0.1800 * x^-0.2

# Exact target function h(x) = (x^2 + 2x)^0.9
exacts <- function(x) (x^2 + 2*x)^0.9

# Randomly sample 100 x-values from the safe domain [40, 240]
# Randomization avoids systematic x‑selection and supports paired inference
set.seed(789)
x <- seq(40, 240, 1)
x100 <- sample(x, size = 100, replace = FALSE)
#
# Sampling nearly half of the available evaluation points while still
# spanning the full 40–240 domain. Using replace = FALSE avoids local
# clustering artifacts, preserves the global geometry of the bias curve,
# and maintains enough density to detect monotonic drift, curvature,
# or sheet‑lift behavior in the signed bias.

# Evaluate GBE2 approximation and exact values at the sampled x
gbe2_approx <- order2(x100)
gbe2_exacts <- exacts(x100)

# signed error
# approx - gbe
signed__error_approx <- (gbe2_approx-gbe2_exacts)
rel_approx_data <- data.frame(x100,signed__error_approx )

ggplot(rel_approx_data,aes(x=x100,y=signed__error_approx)) +
  geom_point(shape = 21, fill = "blue", size = 3) +
  theme_minimal_grid(color="black") +
  ggtitle("Signed error of GBE2 Approximation\nSample size = 100\nGBE2 underestimates function") +
  xlab("x") +
  ylab("Signed Error = GBE approx - Exact") +
  scale_x_continuous(limits=c(35,245),
                              breaks=seq(40,240,20)) +
   geom_hline(yintercept = 0, color = "gray50",linewidth=2.25)


# Signed bias preserves direction (under‑ vs over‑estimation)
gbe2_biases <- gbe2_approx - gbe2_exacts

# gbe2 bias analysis data

gbe2_analysis_data <- data.frame(x100,gbe2_biases)

# ridgeline plot

ggplot(gbe2_analysis_data , aes(x = gbe2_biases, y = "GBE2 Bias")) +
  geom_density_ridges(alpha = 0.6,
                      fill="green") +
  theme_cowplot() +
  scale_x_continuous(limits=c(-0.0025,0.0005)) +
  ggtitle("GBE2 Bias Distribution\nSample size = 100\nMedian = - 0.00031\nGBE2 underestimates ") +
  ylab("Ridge height = Local Density") +
  xlab("Bias") +
  geom_vline(
    xintercept = c(median(gbe2_analysis_data$gbe2_biases),0),
    color = "black",
    linetype = "dashed",
    linewidth =1.25   )

# Evidence of skew
library(moments)
agostino.test(gbe2_analysis_data$gbe2_biases)

# Direction of skewness
library(e1071)
skewness(gbe2_analysis_data$gbe2_biases)


# plot gbe2_biases vs x100
ggplot(gbe2_analysis_data, aes(x = x100, y = gbe2_biases)) +
  geom_point(alpha = 1, color = "blue",size=2.5) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black",
    fill = "steelblue",   # SE ribbon color
    linewidth = 1.2,
    alpha=0.6
  ) +

  theme_cowplot() +
  xlab("x") +
  ylab("GBE2 Bias (approx - exact)") +
  ggtitle("GBE2 Bias vs x with Linear Trend\nRegression quantifies the bias correction") +
  scale_x_continuous(limits=c(35,245),
                     breaks=seq(40,240,20)) 

# lm gbe2_biases vs x100

gbe2_lm <-lm( gbe2_analysis_data$gbe2_biases~gbe2_analysis_data$x100)
summary(gbe2_lm)

# slope and intercepts
coef(gbe2_lm)

# gbe2_corr_function
gbe2_corr <- function(x) x^1.8 +
1.8000 * x^0.8 -  0.1800 * x^-0.2 -  (4.907311e-06*x -1.175748e-03  )

# gbe2_corr values
gbe2_corr_values <- gbe2_corr(x100)

# gbe2_corr biases
gbe2_corr_biases <- gbe2_corr_values - exacts(x100)

# add gbe2_corr biases to gbe2_analysis_data
gbe2_analysis_data$gbe2_corr_biases <-gbe2_corr_biases

# Convert gbe2_analysis_data format to long

library(dplyr)
library(tidyr)

bias_long <- gbe2_analysis_data %>%
  select(all_of(c("gbe2_biases", "gbe2_corr_biases"))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "type",
    values_to = "bias"  )

library(dplyr)
library(tidyr)

median_df <- bias_long %>%
  group_by(type) %>%
  summarize(median_bias = median(bias), .groups = "drop")


library(ggridges)
library(ggplot2)
library(cowplot)

ggplot(bias_long, aes(x = bias, y = type, fill = type)) +
  geom_density_ridges(alpha = 0.6) +
  geom_vline(
    data = median_df,
    aes(xintercept = median_bias, color = type),
    linewidth = 1.2,
    linetype = "dashed",
    show.legend = FALSE
  ) +
  theme_cowplot() +
  scale_fill_manual(values = c("darkgreen", "steelblue")) +
  scale_color_manual(values = c("darkgreen", "steelblue")) +
  xlab("Signed Error") +
  ylab("Ridge height = Local Density") +
  ggtitle("Raw vs Corrected GBE2 Signed Error\nMedian shift reveals collapse toward zero")


# paired wilcoxon test
wilcox.test(
  gbe2_analysis_data$gbe2_biases,
  gbe2_analysis_data$gbe2_corr_biases,
  paired = TRUE,
  exact = FALSE
)

# box plot generation

box_1 <- ggplot(bias_long, aes(x = bias, y = type, fill = type)) +
  geom_boxplot(notch = TRUE) +
  theme_cowplot() +
  ggtitle("Box plots of GBE 2 Bias by Type") +
  theme(legend.position = "none") +
  geom_vline(
    xintercept = 0,
    linewidth = 1,
    linetype = "dashed"
  )

box_1 + coord_flip()

# % median change
(0.0000231 - -0.000311 )/-0.000311 

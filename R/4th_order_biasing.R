exact <- function(x) (x^2 + 2*x)^0.9
approx <- function(x) x^1.8 +
  1.8000 * x^0.8 -
  0.1800 * x^-0.2 +
  0.1320 * x^-1.2 -
  0.1386 * x^-2.2

xs <- seq(40,240,2)



diffs <- approx(xs) - exact(xs)

sample_data <- data.frame(xs,diffs,approx(xs),exact(xs))

# get random data of size 100 from sample data
set.seed(123)
idx <- sample(seq_len(nrow(sample_data)), size = 100, replace = FALSE)
sample_diffs <- sample_data[idx, ]

# perform linear reression on diffs and xs
fit <- lm(sample_diffs$diffs~sample_diffs$xs)
summary(fit)

rel <- diffs / exact(xs)
summary(lm(rel ~ xs))

library(ggplot2)
library(cowplot)

plot <- ggplot(sample_diffs,aes(x=xs,y=diffs)) +
  geom_point(shape = 21, fill = "black", size = 3) 

plot + stat_smooth(method=lm,fill="blue",alpha=0.3,
                   color="black")+theme_minimal_grid(font_size = 16,color="black")+
  ggtitle("Bias vs x\nBias = Approx - True Value")+
  scale_x_continuous(limits=c(30,250),breaks=seq(30,240,30)) +
  ylab("Bias") +
  xlab("x") +
  
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")) +
  annotate("text", x = 55, 
           y = -1.2e-6, label = "- 1.24 e -06",fontface="bold",
           size=5) +
  annotate("text", x = 240, 
         y = -7e-08, label = "- 3.12 e - 09",fontface="bold",
         size=5)



ggplot(sample_diffs,aes(x=diffs)) +
  geom_boxplot(notch = TRUE, fill = "green",outlier.size = 4) +
  theme_cowplot() +
  xlab("Bias") +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")) +
  ggtitle("Bias Box Plot\nBias = GBE Approx - True Value\nSample Size = 101\n15 Outliers")


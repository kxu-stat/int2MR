# Load required libraries
library(ggplot2)  # for plotting
library(dplyr)    # for data manipulation

# Filter results for FDR < 0.1 and only the 3-sample analyses, 
# then ensure 'exposure' is a character vector
res <- load("~/int2MR_ADHD_result/ADHD_result.csv")
res <- res %>%
  filter(FDR_adjusted_pval_beta_int < 0.1) %>%
  filter(samples == "3") %>%
  mutate(exposure = as.character(exposure))

# Prepare male estimates: select exposure, estimate, and SE, 
# then tag as "Male"
res_reshape1 <- res %>%
  select(exposure, est_beta, se_beta) %>%
  mutate(Sex = "Male")

# Prepare female estimates: select exposure, combined estimate (beta + int) and SE, 
# then tag as "Female"
res_reshape2 <- res %>%
  select(exposure, est_beta_plus_int, se_beta_plus_int) %>%
  mutate(Sex = "Female")

# Rename female columns to match the male version for binding
colnames(res_reshape2) <- colnames(res_reshape1)

# Combine male and female data into one dataset
res_reshape <- bind_rows(res_reshape1, res_reshape2)

# Compute 95% confidence interval bounds
res_reshape <- res_reshape %>%
  mutate(
    lower = est_beta - 1.96 * se_beta,  # lower CI
    upper = est_beta + 1.96 * se_beta   # upper CI
  )

# Create the plot: point estimates with error bars, dodge by sex
plot_res <- ggplot(res_reshape, aes(x = est_beta, y = exposure, color = Sex)) +
  geom_point(size = 3, position = position_dodge(0.5)) +                   # points
  geom_errorbar(aes(xmin = lower, xmax = upper),                           # error bars
                position = position_dodge(0.5),
                size = 1, width = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1) +              # vertical zero line
  scale_x_continuous(limits = c(-1, 1), expand = c(0, 0)) +               # x-axis limits
  labs(
    x       = NULL,                                                       # no x-axis label
    y       = NULL,                                                       # no y-axis label
    caption = "Effect Size (+/- 95% CI)"                                   # caption
  ) +
  theme_minimal(base_size = 16) +                                          # minimal theme
  theme(
    legend.position = "right",                                             # legend on the right
    legend.title    = element_text(size = 14),                             # legend title size
    legend.text     = element_text(size = 12),                             # legend text size
    axis.text.y     = element_text(size = 12),                             # y-axis text size
    plot.title      = element_text(hjust = 0.5),                           # center title
    plot.caption    = element_text(hjust = 0.5, size = 16, margin = margin(t = 10))  # caption style
  )

# Display the plot
plot_res

# ==============================================================================
# Script: Plot Expected Null Distribution as Barplot with 95% CI
# ==============================================================================
# Ensure ggplot2 and scales are installed: 
# install.packages(c("ggplot2", "scales"))

library(ggplot2)
library(scales)

# 1. Load the simulated null distribution data
input_file <- "expected_null_distribution_GTR_Gamma.csv"
dist_df <- read.csv(input_file)

# 2. Data Cleaning for Log10 Scale
# Barplots originating from 0 will fail on a log10 scale.
# We set a small pseudo-count (e.g., 0.1) for values <= 0 to ensure proper rendering.
dist_df$mean[dist_df$mean <= 0] <- 0.1
dist_df$lower_95[dist_df$lower_95 <= 0] <- 0.1
dist_df$upper_95[dist_df$upper_95 <= 0] <- 0.1

# 3. Create the Publication-Ready Barplot
# Note: x is converted to factor so bars are drawn as discrete categories
p <- ggplot(dist_df, aes(x = as.factor(simulated_homoplasy_count), y = mean)) +
  
  # Add the bars representing the mean expected value
  geom_bar(stat = "identity", fill = "#3498db", color = "black", 
           width = 0.7, alpha = 0.85, linewidth = 0.5) +
  
  # Add the 95% Confidence Interval error bars
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), 
                width = 0.25, linewidth = 0.8, color = "black") +
  
  # Set Y-axis to Log10 scale for better visualization of rare homoplasy events
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  
  # Add academic labels
  labs(
    x = "Homoplasy Count (Independent Convergent Events)",
    y = "Number of Sites",
    title = "Expected Null Distribution of Convergent Mutations",
    subtitle = "Bars represent mean expected sites; Error bars indicate 95% CI (1,000 Replicates)"
  ) +
  
  # Apply a clean, classic theme suitable for Lancet-style journals
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey30"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    # Remove legend since the subtitle explains the aesthetic
    legend.position = "none" 
  )

# 4. Save the plot as a high-resolution PDF and PNG
ggsave("Figure_Null_Distribution_Barplot.pdf", plot = p, width = 8, height = 6, dpi = 300)
ggsave("Figure_Null_Distribution_Barplot.png", plot = p, width = 8, height = 6, dpi = 300)

cat("Barplot successfully generated and saved as PDF and PNG!\n")
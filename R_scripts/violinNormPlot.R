# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyverse)
library(stringr)
library(ggsignif)  # For adding statistical significance
library(multcompView)
library(car)
library(effectsize)


# Load dataset
genetable_normdata <- read_csv("../data/genetable_normdata.csv")


# Filter the data set by AMR
AMRTable_filtered <- genetable_normdata %>%
  filter(grepl("card", DATABASE)) %>%
  #filter(grepl("vfdb", DATABASE)) %>%
  #filter(grepl("plasmidfinder", DATABASE)) %>%
  mutate(TPM_log_count = log(TPM + 1))  # Adding 1 to avoid log(0)

# Aggregate TPM_log_count values per sample    
AMRTable_aggregated <- AMRTable_filtered %>%
  group_by(sample, Treatment) %>%
  summarise(mean_TPM_log = mean(TPM_log_count, na.rm = TRUE), .groups = "drop") 
    

# one-way ANOVA analysis
anova_result <- aov(mean_TPM_log ~ Treatment, data = AMRTable_aggregated)
summary(anova_result)

# Calculate eta squared
eta_squared(anova_result, partial = FALSE)   # η² (classical)
eta_squared(anova_result, partial = TRUE)    # partial η²

# η²: proportion of total variance explained by Treatment.
# Partial η²: proportion of variance explained by Treatment relative to error + Treatment.

# η² ≈ 0.28–0.29 (Treatment explains ~28–29% of variance in ARG abundance).

# This is considered a large effect size in ecological/biological contexts.

# Extract residuals from the ANOVA model
residuals_anova <- residuals(anova_result)

# Shapiro-Wilk test for normality
shapiro.test(residuals_anova)

# If p-value > 0.05 → residuals are approximately normal.
# If p-value < 0.05 → consider non-parametric alternatives (e.g., Kruskal–Wallis).

# Interpretation
# The null hypothesis of the Shapiro–Wilk test is that the residuals are normally distributed.

# Since p > 0.05, you fail to reject the null hypothesis.

# This means your ANOVA residuals are consistent with normality — the assumption of normality is satisfied.

# Check homogeneity of variance with Levene’s test
leveneTest(mean_TPM_log ~ Treatment, data = AMRTable_aggregated)

# If p > 0.05 → variances are homogeneous, ANOVA assumptions are met.

# If p < 0.05 → variances differ significantly, and you should consider Welch’s ANOVA or a non-parametric Kruskal–Wallis test.

# Levene’s test Interpretation
# The null hypothesis of Levene’s test is that group variances are equal.

# Since p > 0.05, you fail to reject the null hypothesis.

# This means your data does not show evidence of unequal variances — the assumption of homogeneity of variance is satisfied.

# Summary of assumption checks
# Normality (Shapiro–Wilk): p = 0.08196 → residuals are approximately normal.

# Homogeneity of variance (Levene’s): p = 0.1043 → variances are homogeneous.

# 👉 Both assumptions are met, so your one-way ANOVA is valid for comparing AMR profiles across treatments.


# Kruskal–Wallis test (non-parametric alternative to one-way ANOVA)
kruskal_result <- kruskal.test(mean_TPM_log ~ Treatment, data = AMRTable_aggregated)
print(kruskal_result)
kruskal_result$p.value
# Passing to ggplot
label_KW <- paste("Kruskal-Wallis p =", format(kruskal_result$p.value, digits = 2))

# Extract values
H <- kruskal_result$statistic
k <- length(unique(AMRTable_filtered$Treatment))
# Correct total sample size: number of unique samples
n <- length(unique(AMRTable_filtered$sample))

# Compute eta squared for Kruskal–Wallis
eta2_kw <- (H - k + 1) / (n - k)

cat("\n--- Kruskal–Wallis Effect Size (η²_KW) ---\n")
print(eta2_kw)

# Add η² to your ggplot annotation (optional)
label_eta2 <- paste("η² (KW) =", round(eta2_kw, 3))

# Define all 3 pairwise comparisons
my_comparisons <- list( 
  c("Reference diet", "Dulce"),
  c("Soyabean meal", "Reference diet"),
  c("Dulce", "Soyabean meal")
  
)


# Generate the violin plot with statistical significance for AMR
ggplot(AMRTable_filtered, aes(x = Treatment, y = TPM_log_count, fill = Treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  # Jittered points for individual samples
  #geom_jitter(width = 0.1, size = 2, alpha = 0.3) +
  scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                               "Dulce" = "#2ca02c",
                               "Soyabean meal" = "#ff7f0e")) +
  labs(x = "Sample groups", y = "ARGs abundance (log(TPM))") +
  theme_classic() +
  theme(legend.position = "none", text = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust = 0.5)
  ) +
  # Add the pairwise comparisons using Wilcoxon test
  stat_compare_means(comparisons = my_comparisons, 
                     label.y = c(
                       max(AMRTable_filtered$TPM_log_count) + 0.5,  # Control_4 vs Control_5
                       max(AMRTable_filtered$TPM_log_count) + 1.3,  # Control_4 vs Seaweed_4
                       max(AMRTable_filtered$TPM_log_count) + 1.8  # Control_5 vs Seaweed_5
                     ), 
                     method = "wilcox.test",
                     label = "p.signif") + # Shows asterisks (*, **, ***)
    #geom_text(aes(label = Letters), vjust = -0.5, size = 5) +
  annotate("text", x = 2, y = max(AMRTable_filtered$TPM_log_count) + 3,
           #label = paste("ANOVA p =", format(summary(anova_result)[[1]][["Pr(>F)"]][1], digits = 2)),
           label = label_KW,
           size = 4, color = "black") 
  #coord_flip()
  
  
  
  
  # CLEAN UP #################################################
  
  # Clear environment
  rm(list = ls()) 
  
  # Clear plots
  if (length(dev.list()) > 0) {
    dev.off()  # Only if there is an active plot
  }
  
  # Clear console
  cat("\014")  # ctrl+L
  

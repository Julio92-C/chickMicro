# Purpose: Violin + bar plots of species-level read counts across treatment groups
#          from the abricate+Bracken merged dataset, with Kruskal-Wallis test.
# Input:   data/abri_kraken2Bracken_merged.csv (cols: sample, name, sampleCount,
#                                               Treatment, ...)
# Output:  Plots rendered to the active graphics device.

library(readr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(stringr)

# Load dataset
abri_kraken2_clean <- read_csv("data/abri_kraken2Bracken_merged.csv")

# Filter to species-level taxa with sufficient counts
class_data <- abri_kraken2_clean %>%
  filter(!grepl("Terrabacteria group|Bacteroidota/Chlorobiota group|FCB group", name)) %>%
  filter(!(str_count(name, "\\S+") == 1 & name != "Enterococcus")) %>%
  filter(sampleCount > 10) %>%
  rename(count = sampleCount) %>%
  mutate(log_count = log(count + 1))

# Summary statistics
summary_stats <- class_data %>%
  group_by(Treatment) %>%
  summarise(mean_count = mean(count),
            se_count   = sd(count) / sqrt(n()))

classTable_aggregated <- class_data %>%
  group_by(sample, Treatment) %>%
  summarise(mean_log_count = mean(log_count, na.rm = TRUE), .groups = "drop")

# Kruskal-Wallis test
kruskal_result <- kruskal.test(mean_log_count ~ Treatment, data = classTable_aggregated)
print(kruskal_result)
label_KW <- paste("Kruskal-Wallis p =", format(kruskal_result$p.value, digits = 2))

# Per-sample error bars
sample_stats <- class_data %>%
  group_by(sample) %>%
  summarise(mean_count = mean(count),
            se_count   = sd(count) / sqrt(n()),
            Treatment  = first(Treatment))

# Bar plot per sample
Species_count <- ggplot(class_data, aes(x = sample, y = count, fill = Treatment)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(count / 1, 1)), vjust = -0.6, size = 3) +
  geom_hline(aes(yintercept = mean(count)), linetype = "dashed", color = "red") +
  annotate("text", x = Inf, y = mean(class_data$count), label = "",
           hjust = 1.1, vjust = -0.5, color = "red", size = 3) +
  geom_line(aes(group = 1), color = "black", linewidth = 1) +
  geom_errorbar(data = sample_stats,
                aes(x = sample, ymin = mean_count - se_count, ymax = mean_count + se_count),
                width = 0.2, color = "white", inherit.aes = FALSE) +
  geom_point(data = sample_stats, aes(x = sample, y = mean_count),
             color = "red", size = 2, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                               "Dulce"          = "#2ca02c",
                               "Soyabean meal"  = "#ff7f0e")) +
  labs(
    title = bquote(italic("Enterobacteriaceae") ~ "Count per Treatment (Kruskal-Wallis p =" ~
                     .(signif(kruskal_result$p.value, 3)) ~ ")"),
    x = "Sample groups",
    y = "Count"
  ) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        plot.title  = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Treatment, scales = "free_x", nrow = 1)

Species_count

# Violin plot with pairwise comparisons
my_comparisons <- list(
  c("Reference diet", "Dulce"),
  c("Soyabean meal",  "Reference diet"),
  c("Dulce",          "Soyabean meal")
)

ggplot(class_data, aes(x = Treatment, y = log_count, fill = Treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                               "Dulce"          = "#2ca02c",
                               "Soyabean meal"  = "#ff7f0e")) +
  labs(x = "Sample groups", y = "Species abundance (log10(Count))") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(max(class_data$log_count) + 0.5,
                                 max(class_data$log_count) + 1.3,
                                 max(class_data$log_count) + 1.8),
                     method = "wilcox.test",
                     label  = "p.signif") +
  annotate("text", x = 2, y = max(class_data$log_count) + 3,
           label = label_KW, size = 4, color = "black")

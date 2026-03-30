# Purpose: Violin and bar plots of Shannon diversity (or Richness) across
#          treatment groups from the wf-metagenomics diversity output, with
#          Kruskal-Wallis test annotation.
# Input:   data/wf-metagenomics-diversity.csv (rows: Indices, cols: samples)
#          data/chicken_metadata.csv          (cols: sample, Treatment, ...)
# Output:  Plots rendered to the active graphics device.

library(readr)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(plotly)
library(paletteer)
library(funrar)

# Load datasets
wf_metagenomics_diversity <- read_csv("data/wf-metagenomics-diversity.csv")
chicken_metadata <- read_csv("data/chicken_metadata.csv")
chicken_metadata <- na.omit(chicken_metadata)

# Standardise treatment labels
modify_strings <- function(df, column_name, target_string, replacement_string) {
  df %>% mutate(!!sym(column_name) :=
                  str_replace_all(!!sym(column_name), fixed(target_string), replacement_string))
}

chicken_metadata <- modify_strings(chicken_metadata, "Treatment",
                                   "Reference diet (Wheat -soyabean)", "Reference diet")
chicken_metadata <- modify_strings(chicken_metadata, "Treatment",
                                   "Seaweed", "Dulce")

# Pivot the diversity dataset
diversity_pivoted <- wf_metagenomics_diversity %>%
  pivot_longer(cols = -Indices, names_to = "sample", values_to = "value") %>%
  pivot_wider(names_from = Indices, values_from = value)

# Remove last row (Total count row)
diversity_pivoted <- diversity_pivoted[-19, ]

# Merge with metadata
diversity_merged <- merge(diversity_pivoted, chicken_metadata, by = "sample")

# Select the index of interest (swap comments to switch to Richness)
diversity_filtered <- diversity_merged %>%
  arrange(Treatment) %>%
  select(sample, Treatment, `Shannon diversity index`) %>%
  rename(Index = `Shannon diversity index`)

# Summary statistics
summary_stats <- diversity_filtered %>%
  group_by(Treatment) %>%
  summarise(mean_count = mean(Index),
            se_count   = sd(Index) / sqrt(n()))

# Aggregate per sample
classTable_aggregated <- diversity_filtered %>%
  group_by(sample, Treatment) %>%
  summarise(mean_Indece = mean(Index, na.rm = TRUE), .groups = "drop")

# Kruskal-Wallis test
kruskal_result <- kruskal.test(mean_Indece ~ Treatment, data = classTable_aggregated)
print(kruskal_result)
label_KW <- paste("Kruskal-Wallis p =", format(kruskal_result$p.value, digits = 2))

# Pairwise comparisons
my_comparisons <- list(
  c("Reference diet", "Dulce"),
  c("Soyabean meal",  "Reference diet"),
  c("Dulce",          "Soyabean meal")
)

# Violin plot
ggplot(diversity_filtered, aes(x = Treatment, y = Index, fill = Treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                               "Dulce"          = "#2ca02c",
                               "Soyabean meal"  = "#ff7f0e")) +
  labs(x = "Sample groups", y = "Richness") +
  theme_classic() +
  theme(legend.position = "top",
        text = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  annotate("text", x = 2, y = max(diversity_filtered$Index) + 30,
           label = label_KW, size = 5, color = "black")

# Per-sample mean for bar plot
sample_stats <- diversity_filtered %>%
  group_by(sample) %>%
  summarise(mean_count = mean(Index),
            se_count   = sd(Index) / sqrt(n()),
            Treatment  = first(Treatment))

# Bar plot
Species_count <- ggplot(diversity_filtered, aes(x = sample, y = Index, fill = Treatment)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(Index, 1)), vjust = -0.6, size = 3) +
  geom_hline(aes(yintercept = mean(Index)), linetype = "dashed", color = "red") +
  annotate("text", x = Inf, y = mean(diversity_filtered$Index), label = "",
           hjust = 1.1, vjust = -0.5, color = "red", size = 3) +
  geom_line(aes(group = 1), color = "black", linewidth = 1) +
  geom_point(data = sample_stats, aes(x = sample, y = mean_count),
             color = "red", size = 2, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                               "Dulce"          = "#2ca02c",
                               "Soyabean meal"  = "#ff7f0e")) +
  labs(
    title = bquote(italic("") ~ "Shannon index by Treatment (Kruskal-Wallis p =" ~
                     .(signif(kruskal_result$p.value, 3)) ~ ")"),
    x = "Sample groups",
    y = "Diversity"
  ) +
  theme_classic() +
  theme(legend.position = "none",
        text            = element_text(size = 14),
        axis.text.y     = element_text(angle = 90, hjust = 0.5),
        plot.title      = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x     = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Treatment, scales = "free_x", nrow = 1)

Species_count

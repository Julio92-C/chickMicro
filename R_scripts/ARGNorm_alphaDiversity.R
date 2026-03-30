# Purpose: Alpha-diversity (Richness and Shannon) of MGE/ARG/VF gene profiles
#          across treatment groups, with Kruskal-Wallis significance annotation.
# Input:   data/genetable_normdata.csv (cols: sample, GENE, DATABASE, TPM, Treatment, ...)
# Output:  Violin plot and bar plot rendered to the active graphics device.

library(vegan)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(plotly)
library(paletteer)
library(funrar)

# Load dataset
genetable_normdata <- read_csv("data/genetable_normdata.csv")

# Metadata
metadata <- genetable_normdata %>%
  distinct(sample, .keep_all = TRUE) %>%
  select(sample, Treatment)

# Filter to one element type (uncomment the desired line)
geneTable_filtered <- genetable_normdata %>%
  #filter(grepl("card", DATABASE)) %>%
  #filter(grepl("vfdb", DATABASE)) %>%
  filter(grepl("plasmidfinder", DATABASE)) %>%
  mutate(TPM_log_count = log(TPM + 1)) %>%
  select(sample, GENE, TPM_log_count) %>%
  tidyr::pivot_wider(
    names_from  = GENE,
    values_from = TPM_log_count,
    values_fill = 0,
    values_fn   = sum
  ) %>%
  distinct(sample, .keep_all = TRUE) %>%
  column_to_rownames("sample")

abundance_matrix <- as.matrix(geneTable_filtered)

# Alpha-diversity metrics
arg_richness <- specnumber(abundance_matrix)
arg_shannon  <- diversity(abundance_matrix, index = "shannon")

alpha_results <- data.frame(
  sample  = rownames(abundance_matrix),
  Richness = as.numeric(arg_richness),
  Shannon  = as.numeric(arg_shannon)
)

alphaMeta_results <- merge(alpha_results, metadata, by = "sample")

# Select index to plot (switch between Richness / Shannon here)
diversity_filtered <- alphaMeta_results %>%
  arrange(Treatment) %>%
  select(sample, Treatment, Richness) %>%
  rename(Index = Richness)
  # select(sample, Treatment, Shannon) %>%
  # rename(Index = Shannon)

# Summary statistics
summary_stats <- diversity_filtered %>%
  group_by(Treatment) %>%
  summarise(mean_Index = mean(Index),
            se_Index   = sd(Index) / sqrt(n()))

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
  theme(legend.position = "none",
        text = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  annotate("text", x = 2, y = max(diversity_filtered$Index) + 5,
           label = label_KW, size = 5, color = "black")

# Per-sample bar plot
sample_stats <- diversity_filtered %>%
  group_by(sample) %>%
  summarise(mean_count = mean(Index),
            se_count   = sd(Index) / sqrt(n()),
            Treatment  = first(Treatment))

Species_count <- ggplot(diversity_filtered, aes(x = sample, y = Index, fill = Treatment)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(Index / 1, 1)), vjust = -0.6, size = 3) +
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
    title = bquote(italic("") ~ "Richness by Treatment (Kruskal-Wallis p =" ~
                     .(signif(kruskal_result$p.value, 3)) ~ ")"),
    x = "Sample groups",
    y = "Richness"
  ) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        plot.title  = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Treatment, scales = "free_x", nrow = 1)

Species_count

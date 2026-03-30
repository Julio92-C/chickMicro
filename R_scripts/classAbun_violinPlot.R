# Purpose: Violin plot and Kruskal-Wallis test of log-transformed class-level
#          read counts (Clostridia, Bacilli, etc.) across treatment groups
#          from the Bracken output.
# Input:   data/bracken_arranged.csv (cols: #perc, tot_all, tot_lvl, taxid, name, samples...)
#          data/chicken_metadata.csv (cols: sample, Treatment, ...)
# Output:  Violin plot rendered to the active graphics device.

library(ggplot2)
library(dplyr)
library(readr)
library(tidyverse)
library(stringr)
library(ggsignif)
library(multcompView)
library(car)
library(effectsize)

# Load datasets
bracken_arranged <- read_csv("data/bracken_arranged.csv")
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

# Pivot longer and filter to target classes
bracken_pivoted <- bracken_arranged %>%
  pivot_longer(cols = 4:21, names_to = "sample", values_to = "count") %>%
  filter(!grepl("root|Homo sapiens|cellular organisms|unclassified|Bacteria|environmental samples", name)) %>%
  filter(grepl("Clostridia|Actinomycetes|Bacilli|Gammaproteobacteria|Alphaproteobacteria|Betaproteobacteria|Bacteroidia|Spirochaetia", name)) %>%
  filter(name != "Clostridiaceae") %>%
  filter(!is.na(name)) %>%
  group_by(sample) %>%
  filter(count > 5) %>%
  mutate(Percentage = (count / sum(count)) * 100) %>%
  ungroup()

# Merge with metadata
bracken_merged <- merge(bracken_pivoted, chicken_metadata, by = "sample")

uniqueSpecies <- length(unique(bracken_merged$name))

# Log transform counts
classTable_filtered <- bracken_merged %>%
  mutate(log_count = log(count + 1))

# Summary statistics
summary_stats <- classTable_filtered %>%
  group_by(Treatment) %>%
  summarise(mean_count = mean(count),
            se_count   = sd(count) / sqrt(n()))

# Aggregate per sample
classTable_aggregated <- classTable_filtered %>%
  group_by(sample, Treatment) %>%
  summarise(mean_log_count = mean(log_count, na.rm = TRUE), .groups = "drop")

# Kruskal-Wallis test
kruskal_result <- kruskal.test(mean_log_count ~ Treatment, data = classTable_aggregated)
print(kruskal_result)
label_KW <- paste("Kruskal-Wallis p =", format(kruskal_result$p.value, digits = 2))

# Kruskal-Wallis effect size (eta squared)
H        <- kruskal_result$statistic
k        <- length(unique(classTable_aggregated$Treatment))
n        <- length(unique(classTable_aggregated$sample))
eta2_kw  <- (H - k + 1) / (n - k)
cat("\n--- Kruskal-Wallis Effect Size (eta2_KW) ---\n")
print(eta2_kw)
label_eta2 <- paste("eta2 (KW) =", round(eta2_kw, 3))

# Pairwise comparisons
my_comparisons <- list(
  c("Reference diet", "Dulce"),
  c("Soyabean meal",  "Reference diet"),
  c("Dulce",          "Soyabean meal")
)

# Violin plot
ggplot(classTable_filtered, aes(x = Treatment, y = log_count, fill = Treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                               "Dulce"          = "#2ca02c",
                               "Soyabean meal"  = "#ff7f0e")) +
  labs(x = "Sample groups", y = "Class abundance (log10(Count))") +
  theme_classic() +
  theme(legend.position = "none",
        text        = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  annotate("text", x = 2, y = max(classTable_filtered$log_count) + 3.0,
           label = label_KW, size = 4, color = "black")

# Purpose: Stacked relative-abundance bar plot of top-50 taxa detected across
#          samples in the abricate+Bracken merged dataset.
# Input:   data/abri_kraken2Bracken_merged.csv (cols: sample, name, sampleCount,
#                                               DATABASE, Treatment, ...)
# Output:  Plot rendered to the active graphics device.

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(purrr)
library(paletteer)
library(plotly)
library(htmlwidgets)
library(stringr)

# Load dataset
abri_kraken2_merged <- read_csv("data/abri_kraken2Bracken_merged.csv")

# Filter and abbreviate species names
abri_kraken2_filtered <- abri_kraken2_merged %>%
  filter(grepl("card|vfdb|plasmidfinder", DATABASE)) %>%
  filter(!grepl("Terrabacteria group|Bacteroidota/Chlorobiota group|FCB group", name)) %>%
  filter(str_count(name, "\\S+") > 1) %>%
  mutate(
    name = str_replace(name, "^(\\w)\\w+\\s(\\w+)", "\\1. \\2"),
    name = if_else(
      str_length(name) > 30,
      str_replace(name, "^((?:\\S+\\s+){2}).*", "\\1"),
      name
    )
  ) %>%
  mutate(
    name = case_when(
      grepl("\\[Clostridium\\] asparagiforme ", name) ~ "C. asparagiforme",
      grepl("\\[Clostridium\\] scindens", name)       ~ "C. scindens",
      grepl("\\[Ruminococcus\\] lactaris", name)      ~ "R. lactaris ",
      grepl("\\[Ruminococcus\\] torques", name)       ~ "R. torques",
      grepl("\\[Clostridium] innocuum", name)         ~ "C. innocuum",
      TRUE ~ name
    ),
    name = ifelse(grepl("u. Subdoligranulum", name), "U. Subdoligranulum", name)
  ) %>%
  mutate(Percentage = (sampleCount / sum(sampleCount)) * 100) %>%
  ungroup()

# Identify top-50 taxa by total abundance
top50_taxa <- abri_kraken2_filtered %>%
  group_by(name) %>%
  summarise(total = sum(sampleCount, na.rm = TRUE)) %>%
  arrange(desc(total)) %>%
  slice_head(n = 50) %>%
  pull(name)

# Collapse remaining taxa into "Others"
abri_kraken2_collapsed <- abri_kraken2_filtered %>%
  mutate(name = ifelse(name %in% top50_taxa, name, "Others")) %>%
  group_by(sample) %>%
  mutate(PercTop50 = (sampleCount / sum(sampleCount)) * 100) %>%
  arrange(sample, desc(PercTop50)) %>%
  ungroup()

abri_kraken2_collapsed$name <- factor(
  abri_kraken2_collapsed$name,
  levels = abri_kraken2_collapsed %>%
    group_by(name) %>%
    summarise(total = sum(PercTop50)) %>%
    arrange(desc(total)) %>%
    pull(name)
)

taxa_f    <- as.factor(abri_kraken2_collapsed$name)
t         <- length(levels(taxa_f))
t_palette <- paletteer_d("ggsci::default_igv", n = t)

ann_name_col <- abri_kraken2_collapsed %>%
  select(name) %>%
  na.omit() %>%
  distinct(name) %>%
  mutate(color = t_palette)

name_colors <- setNames(as.character(ann_name_col$color), ann_name_col$name)

# Stacked relative-abundance bar plot
ra <- ggplot(abri_kraken2_collapsed,
             aes(x = factor(sample), y = PercTop50, fill = name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = t_palette) +
  labs(x = "Sample Groups by Treatment", y = "Percentage (Relative abundance)",
       fill = "Taxonomic level") +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14)) +
  guides(fill = guide_legend(nrow = 10)) +
  facet_wrap(~ Treatment, scales = "free_x", nrow = 1)

ra

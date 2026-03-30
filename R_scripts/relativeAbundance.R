# Purpose: Stacked relative-abundance bar plot of taxa detected across samples,
#          faceted by treatment group (Reference diet, Seaweed, Soyabean meal).
# Input:   ../data/abri_kraken2_filtered.csv (cols: SAMPLE, NAME, TAXID, GENE,
#                                             DATABASE, COVERAGE_PCT, IDENTITY_PCT, ...)
#          ../data/chicken_metadata.csv       (cols: SAMPLE, TREATMENT)
# Output:  Stacked bar chart rendered to the active graphics device.

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

# Load datasets
abri_kraken2_filtered <- read_csv("../data/abri_kraken2_filtered.csv")
chicken_metadata <- read_csv("../data/chicken_metadata.csv")
chicken_metadata <- na.omit(chicken_metadata)

# Merge on SAMPLE
abri_kraken2_merged <- merge(abri_kraken2_filtered, chicken_metadata, by = "SAMPLE")

# Filter and abbreviate species names
abri_kraken2_filtered <- abri_kraken2_merged %>%
  filter(!grepl("Terrabacteria group|Bacteroidota/Chlorobiota group|FCB group", NAME)) %>%
  filter(str_count(NAME, "\\S+") > 1) %>%
  mutate(
    NAME = str_replace(NAME, "^(\\w)\\w+\\s(\\w+)", "\\1. \\2"),
    NAME = if_else(
      str_length(NAME) > 30,
      str_replace(NAME, "^((?:\\S+\\s+){2}).*", "\\1"),
      NAME
    )
  ) %>%
  mutate(
    NAME = case_when(
      grepl("\\[Clostridium\\] asparagiforme ", NAME) ~ "C. asparagiforme",
      grepl("\\[Clostridium\\] scindens", NAME)       ~ "C. scindens",
      grepl("\\[Ruminococcus\\] lactaris", NAME)      ~ "R. lactaris",
      grepl("\\[Ruminococcus\\] torques", NAME)       ~ "R. torques",
      grepl("\\[Clostridium\\] innocuum", NAME)       ~ "C. innocuum",
      TRUE ~ NAME
    ),
    NAME = ifelse(grepl("u. Subdoligranulum", NAME), "U. Subdoligranulum", NAME)
  ) %>%
  group_by(SAMPLE) %>%
  mutate(Percentage = (1 / n()) * 100) %>%
  ungroup()

# Identify top-50 taxa by total count
top50_taxa <- abri_kraken2_filtered %>%
  group_by(NAME) %>%
  summarise(total = n()) %>%
  arrange(desc(total)) %>%
  slice_head(n = 50) %>%
  pull(NAME)

# Collapse remaining taxa into "Others"
abri_kraken2_collapsed <- abri_kraken2_filtered %>%
  mutate(NAME = ifelse(NAME %in% top50_taxa, NAME, "Others")) %>%
  group_by(SAMPLE, NAME, TREATMENT) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(SAMPLE) %>%
  mutate(PercTop50 = (Count / sum(Count)) * 100) %>%
  arrange(SAMPLE, desc(PercTop50)) %>%
  ungroup()

abri_kraken2_collapsed$NAME <- factor(
  abri_kraken2_collapsed$NAME,
  levels = abri_kraken2_collapsed %>%
    group_by(NAME) %>%
    summarise(total = sum(PercTop50)) %>%
    arrange(desc(total)) %>%
    pull(NAME)
)

taxa_f    <- as.factor(abri_kraken2_collapsed$NAME)
t         <- length(levels(taxa_f))
t_palette <- paletteer_d("ggsci::default_igv", n = t)

ann_name_col <- abri_kraken2_collapsed %>%
  select(NAME) %>%
  na.omit() %>%
  distinct(NAME) %>%
  mutate(color = t_palette)

name_colors <- setNames(as.character(ann_name_col$color), ann_name_col$NAME)

# Stacked relative-abundance bar plot
ra <- ggplot(abri_kraken2_collapsed,
             aes(x = factor(SAMPLE), y = PercTop50, fill = NAME)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = t_palette) +
  labs(x = "Sample Groups by Treatment", y = "Percentage (Relative abundance)",
       fill = "Taxonomic level") +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14)) +
  guides(fill = guide_legend(nrow = 10)) +
  facet_wrap(~ TREATMENT, scales = "free_x", nrow = 1)

ra

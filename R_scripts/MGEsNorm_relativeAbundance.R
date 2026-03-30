# Purpose: Stacked relative-abundance bar plots and ranked bar charts for MGE
#          (plasmidfinder replicon) profiles derived from TPM-normalised data.
# Input:   data/abri_kraken2Bracken_merged.csv (cols: sample, name, GENE, PRODUCT,
#                                               DATABASE, Treatment, sequence, ...)
#          data/genetable_normdata.csv          (cols: sample, GENE, DATABASE, TPM,
#                                               Treatment, ...)
# Output:  Plots rendered to the active graphics device.

library(readr)
library(plotly)
library(dplyr)
library(tidyr)
library(stringr)
library(grid)
library(paletteer)
library(ggplot2)
library(ggpubr)

# Load datasets
abri_kraken2Bracken_merged <- read_csv("data/abri_kraken2Bracken_merged.csv")
genetable_normdata          <- read_csv("data/genetable_normdata.csv")

# Replicon classification function
classify_plasmids <- function(df) {
  df %>% mutate(Replicon_Family = case_when(
    str_detect(GENE, "^Col")  ~ "Col-like",
    str_detect(GENE, "^IncF") ~ "IncF",
    str_detect(GENE, "^IncX") ~ "IncX",
    str_detect(GENE, "^Inc")  ~ "Other Inc",
    TRUE                      ~ "Unknown/Other"
  ))
}

# Raw abricate summary (unique sequences)
df_filtered <- abri_kraken2Bracken_merged %>%
  filter(DATABASE == "plasmidfinder")

uniqueSeq  <- length(unique(df_filtered$sequence))
uniqueGene <- length(unique(df_filtered$GENE))

df_filtered   <- classify_plasmids(df_filtered)
df_arrange    <- df_filtered %>%
  distinct(GENE, .keep_all = TRUE) %>%
  group_by(Replicon_Family) %>%
  summarise(count = n())

# TPM-normalised replicon relative abundance
genetable_filtered <- genetable_normdata %>%
  filter(DATABASE == "plasmidfinder") %>%
  classify_plasmids() %>%
  select(sample, GENE, PRODUCT, Replicon_Family, TPM, Treatment) %>%
  group_by(sample) %>%
  mutate(Percentage = (TPM / sum(TPM)) * 100) %>%
  ungroup()

uniqueMGEs     <- length(unique(genetable_filtered$GENE))
uniqueFunction <- length(unique(genetable_filtered$Replicon_Family))
uniqueSamples  <- length(unique(genetable_filtered$sample))

MGEs_f    <- as.factor(genetable_filtered$GENE)
f         <- length(levels(MGEs_f))
f_palette <- paletteer_d("ggsci::default_ucscgb", n = f)

# Stacked relative-abundance bar plot
ra <- ggplot(genetable_filtered, aes(x = factor(sample), y = Percentage, fill = GENE)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = f_palette) +
  labs(x = "Sample Groups by Treatment", y = "Percentage (Relative abundance)") +
  theme_classic() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14)) +
  guides(fill = guide_legend(nrow = 10)) +
  labs(fill = "Plasmid replicons") +
  facet_wrap(~ Treatment, scales = "free_x", nrow = 1)

ra

# Ranked total TPM per replicon
pfdb_filtered <- genetable_filtered %>%
  group_by(GENE) %>%
  summarise(Total_TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total_TPM))

PRs_count <- ggplot(pfdb_filtered, aes(x = reorder(GENE, Total_TPM),
                                        y = Total_TPM, fill = GENE)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = median(pfdb_filtered$Total_TPM), linetype = "dashed", color = "black") +
  geom_text(aes(label = paste0(round(Total_TPM / 1000, 1), "K")),
            vjust = 0.5, hjust = 1.0, size = 3.5) +
  scale_fill_manual(values = f_palette) +
  scale_y_log10() +
  labs(title = "", x = "Plasmid replicons", y = "Total Count") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        text = element_text(size = 14)) +
  coord_flip()

PRs_count

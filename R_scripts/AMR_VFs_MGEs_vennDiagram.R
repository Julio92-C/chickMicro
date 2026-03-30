# Purpose: Venn diagram showing overlap of taxa carrying AMR genes, VFs, and MGEs.
# Input:   data/abri_kraken2Bracken_merged.csv (cols: sample, name, GENE, DATABASE,
#                                               RESISTANCE, Treatment, ...)
# Output:  Venn diagram rendered to the active graphics device.

library(readr)
library(dplyr)
library(tidyverse)
library(paletteer)
library(reshape2)
library(stringr)
library(VennDiagram)

# Load dataset
abri_kraken2_merged <- read_csv("data/abri_kraken2Bracken_merged.csv")

# Filter to species-level taxa with AMR/VF/MGE annotations
abri_kraken2_filtered <- abri_kraken2_merged %>%
  filter(grepl("card|vfdb|plasmidfinder", DATABASE)) %>%
  filter(!grepl("Terrabacteria group|Bacteroidota/Chlorobiota group|FCB group", name)) %>%
  filter(str_count(name, "\\S+") > 1) %>%
  tidyr::pivot_wider(
    names_from  = DATABASE,
    values_from = GENE,
    values_fill = " "
  )

AMR_unique   <- unique(as.factor(abri_kraken2_filtered$card))
MGRs_unique  <- unique(as.factor(abri_kraken2_filtered$plasmidfinder))
Vfs_unique   <- unique(as.factor(abri_kraken2_filtered$vfdb))
uniqueTaxa   <- unique(as.factor(abri_kraken2_filtered$name))

# Build per-taxon summary: which databases are present?
df_wide_ann_rows <- abri_kraken2_filtered %>%
  group_by(name) %>%
  summarise(
    ARGs = paste(unique(card),          collapse = ", "),
    VFs  = paste(unique(vfdb),          collapse = ", "),
    MGEs = paste(unique(plasmidfinder), collapse = ", "),
    .groups = "drop"
  ) %>%
  select(name, ARGs, VFs, MGEs) %>%
  distinct(name, .keep_all = TRUE) %>%
  column_to_rownames("name") %>%
  mutate(across(everything(), ~ if_else(. == " ", NA_character_, .))) %>%
  filter(!if_all(everything(), is.na))

# Convert to binary matrix and create Venn sets
sample2arg_data_matrix <- data.matrix(df_wide_ann_rows, rownames.force = NA) %>%
  replace(is.na(.), 0)

sum(is.na(sample2arg_data_matrix))

sample2arg_data_Bmatrix <- as.matrix((sample2arg_data_matrix > 0) + 0)

sets       <- apply(sample2arg_data_Bmatrix, 2, function(col) which(col == 1))
names(sets) <- colnames(sample2arg_data_Bmatrix)
print(sets)

# Venn diagram
venn.plot <- venn.diagram(
  x              = sets,
  category.names = c("AMR", "VFs", "MGEs"),
  fill           = c(card = "#CD3333FF", vfdb = "#EEAD0EFF", plasmidfinder = "blue"),
  alpha          = 0.5,
  height         = 50,
  width          = 50,
  filename       = NULL,
  cat.cex        = 1.2,
  cex            = 1.8
)
grid.newpage()
grid.draw(venn.plot)

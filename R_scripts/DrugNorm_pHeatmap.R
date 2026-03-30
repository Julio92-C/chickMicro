# Purpose: Pheatmap of TPM-normalised drug-class (CARD) profiles aggregated
#          by drug class per sample, with treatment column annotation.
# Input:   data/genetable_normdata.csv (cols: sample, GENE, RESISTANCE, DATABASE,
#                                       TPM, Treatment, ...)
# Output:  Heatmap rendered to the active graphics device.

library(readr)
library(dplyr)
library(tidyverse)
library(paletteer)
library(reshape2)
library(pheatmap)
library(stringr)

# Load dataset
genetable_normdata <- read_csv("data/genetable_normdata.csv")

# Helper functions
clean_string <- function(x) str_to_title(str_replace_all(x, "_", " "))

classify_resistance <- function(res_string) {
  if (is.null(res_string) || is.na(res_string) || res_string == "") return(NA)
  classes        <- trimws(unlist(strsplit(res_string, ";")))
  classes        <- classes[classes != ""]
  unique_classes <- unique(classes)
  if (length(unique_classes) == 1) return(unique_classes)
  if (length(unique_classes) > 1) return("Multi-drug")
  return(NA)
}

# Build drug-class x sample TPM matrix
geneTable_filtered <- genetable_normdata %>%
  mutate(RESISTANCE = str_replace(
    RESISTANCE,
    "lincosamide;macrolide;streptogramin;streptogramin_A;streptogramin_B", "MLS"
  )) %>%
  mutate(RESISTANCE = sapply(RESISTANCE, clean_string)) %>%
  mutate(RESISTANCE = sapply(RESISTANCE, classify_resistance)) %>%
  mutate(RESISTANCE = str_replace(RESISTANCE, "Mls", "MLS")) %>%
  arrange(Treatment) %>%
  filter(grepl("card", DATABASE)) %>%
  select(sample, RESISTANCE, TPM) %>%
  tidyr::pivot_wider(
    names_from  = sample,
    values_from = TPM,
    values_fill = 0,
    values_fn   = sum
  ) %>%
  distinct(RESISTANCE, .keep_all = TRUE) %>%
  column_to_rownames("RESISTANCE")

geneTable_matrix <- data.matrix(geneTable_filtered, rownames.force = NA)
sum(is.na(geneTable_matrix))

# Sort rows by decreasing row sum
rearrange_matrix <- function(mat) {
  mat[order(rowSums(mat), decreasing = TRUE), ]
}
geneTable_matrix_sorted <- rearrange_matrix(geneTable_matrix)

# Clustering callback
callback <- function(hc, mat) {
  sv   <- svd(t(mat))$v[, 1]
  dend <- reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

# Min-max scale per column
geneTable_scaled <- apply(geneTable_matrix_sorted, 2, function(x) (x - min(x)) / (max(x) - min(x)))
range(geneTable_scaled)
rowSums(geneTable_scaled)

# Column annotation: treatment
geneTable_filtered_transposed <- as.data.frame(t(geneTable_filtered)) %>%
  mutate(sample = rownames(.)) %>%
  select(sample)

Metadata_Treatment <- genetable_normdata %>%
  select(sample, Treatment) %>%
  distinct()

ann_col <- inner_join(geneTable_filtered_transposed, Metadata_Treatment, by = "sample") %>%
  na.omit() %>%
  remove_rownames() %>%
  column_to_rownames(var = "sample")

ann_colors <- list(
  Treatment = c("Reference diet" = "#1f77b4",
                "Dulce"          = "#2ca02c",
                "Soyabean meal"  = "#ff7f0e")
)

# Pheatmap
heatmap_plot <- pheatmap(
  geneTable_scaled,
  display_numbers     = FALSE,
  cluster_cols        = FALSE,
  cluster_rows        = FALSE,
  scale               = "none",
  clustering_callback = callback,
  color               = colorRampPalette(c("#f5f06e", "#3df5e9", "#f53d8d"))(100),
  border_color        = NA,
  show_rownames       = TRUE,
  annotation_col      = ann_col,
  annotation_colors   = ann_colors,
  fontsize_row        = 12,
  fontsize_col        = 12,
  fontsize            = 12
)

heatmap_plot

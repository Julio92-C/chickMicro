# Purpose: Pheatmap of TPM-normalised MGE (plasmidfinder replicon) profiles
#          across samples, with treatment column annotation.
# Input:   data/genetable_normdata.csv (cols: sample, GENE, DATABASE, TPM, Treatment, ...)
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

# Build replicon x sample TPM matrix
geneTable_filtered <- genetable_normdata %>%
  filter(DATABASE == "plasmidfinder") %>%
  arrange(Treatment) %>%
  select(sample, GENE, TPM) %>%
  tidyr::pivot_wider(
    names_from  = sample,
    values_from = TPM,
    values_fill = 0,
    values_fn   = sum
  ) %>%
  distinct(GENE, .keep_all = TRUE) %>%
  column_to_rownames("GENE")

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
  color               = colorRampPalette(c("#ebc3c0", "#bd0606", "#d98c07"))(100),
  border_color        = NA,
  show_rownames       = TRUE,
  annotation_col      = ann_col,
  annotation_colors   = ann_colors,
  fontsize_row        = 12,
  fontsize_col        = 12,
  fontsize            = 12
)

heatmap_plot

# Purpose: Pheatmap of TPM-normalised VF gene profiles (effector delivery system
#          subset) with virulence-function row annotation and treatment column annotation.
# Input:   data/genetable_normdata.csv (cols: sample, GENE, PRODUCT, DATABASE,
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

# Extract virulence function category from PRODUCT description
extract_productFunction <- function(string) {
  match  <- regmatches(string, regexpr("\\- ([^\\(]+) \\(", string))
  result <- gsub("\\- | \\(", "", match)
  return(result)
}

# Build GENE x sample TPM matrix (filter by function if desired)
geneTable_filtered <- genetable_normdata %>%
  filter(grepl("vfdb", DATABASE)) %>%
  mutate(Functions = sapply(PRODUCT, extract_productFunction)) %>%
  arrange(Treatment, Functions) %>%
  # Uncomment ONE of the lines below to restrict to a specific function category:
  #filter(Functions == "Adherence") %>%
  filter(Functions == "Effector delivery system") %>%
  #filter(Functions != "Effector delivery system") %>%
  group_by(sample) %>%
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

# Clustering callback
callback <- function(hc, mat) {
  sv   <- svd(t(mat))$v[, 1]
  dend <- reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

# Min-max scale per column
geneTable_scaled <- apply(geneTable_matrix, 2, function(x) (x - min(x)) / (max(x) - min(x)))
range(geneTable_scaled)
rowSums(geneTable_scaled)

# Column annotation: treatment
geneTable_filtered_transposed <- as.data.frame(t(geneTable_filtered)) %>%
  mutate(sample = rownames(.)) %>%
  select(sample)

Metadata_Treatment <- genetable_normdata %>%
  filter(DATABASE == "vfdb") %>%
  select(sample, Treatment) %>%
  distinct()

ann_col <- inner_join(geneTable_filtered_transposed, Metadata_Treatment, by = "sample") %>%
  na.omit() %>%
  remove_rownames() %>%
  column_to_rownames(var = "sample")

# Row annotation: virulence function category
gene_fun <- genetable_normdata %>%
  filter(DATABASE == "vfdb") %>%
  mutate(Functions = sapply(PRODUCT, extract_productFunction)) %>%
  select(GENE, Functions) %>%
  mutate(Functions = str_replace(
    Functions,
    "Antimicrobial activity/Competitive advantage", "Antimicrobial activity/CA*"
  )) %>%
  distinct() %>%
  column_to_rownames(var = "GENE") %>%
  arrange(Functions)

virulentCategory <- as.factor(gene_fun$Functions)
r                <- length(levels(virulentCategory))
r_palette        <- paletteer_d("ggsci::default_nejm", n = r)

ann_row_funt_col <- gene_fun %>%
  select(Functions) %>%
  na.omit() %>%
  distinct(Functions) %>%
  mutate(color = r_palette)

fun_colors <- setNames(as.character(ann_row_funt_col$color), ann_row_funt_col$Functions)

ann_colors <- list(
  Treatment = c("Reference diet" = "#1f77b4",
                "Dulce"          = "#2ca02c",
                "Soyabean meal"  = "#ff7f0e"),
  Functions = fun_colors
)

# Pheatmap
heatmap_plot <- pheatmap(
  geneTable_scaled,
  display_numbers     = FALSE,
  cluster_cols        = FALSE,
  cluster_rows        = FALSE,
  scale               = "none",
  clustering_callback = callback,
  color               = colorRampPalette(c("#0612bd", "#bbbbbd", "#bd0606"))(100),
  border_color        = NA,
  show_rownames       = TRUE,
  annotation_row      = gene_fun,
  annotation_col      = ann_col,
  annotation_colors   = ann_colors,
  fontsize_row        = 12,
  fontsize_col        = 12,
  fontsize            = 12
)

heatmap_plot

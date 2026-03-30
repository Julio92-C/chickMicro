# Purpose: PCoA of MGE (plasmidfinder) gene profiles using Bray-Curtis distances
#          with Hellinger transformation, PERMANOVA, and PERMDISP tests.
# Input:   data/genetable_normdata.csv (cols: sample, GENE, DATABASE, TPM, type, ...)
# Output:  PCoA scatter plot rendered to the active graphics device.
# Note:    This script uses a 'type' metadata column (not 'Treatment'); update the
#          column name below if your metadata uses a different field.

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

# Filter to MGEs; log-transform TPM; pivot to sample x gene matrix
geneTable_filtered <- genetable_normdata %>%
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

# Hellinger transformation and Bray-Curtis distance
relative_abundance_matrix <- decostand(abundance_matrix, method = "hellinger")

sum(is.na(relative_abundance_matrix))
sum(is.nan(as.matrix(relative_abundance_matrix)))

clean_matrix    <- na.omit(relative_abundance_matrix)
bray_curtis_dist <- vegdist(clean_matrix, method = "bray")

# PCoA
pcoa_result  <- cmdscale(bray_curtis_dist, eig = TRUE, k = 2)
eig_values   <- pcoa_result$eig
var_explained <- eig_values / sum(eig_values) * 100

pcoa_df <- data.frame(
  sample = rownames(clean_matrix),
  PC1    = pcoa_result$points[, 1],
  PC2    = pcoa_result$points[, 2]
)

# Metadata — uses 'type' column; change to 'Treatment' if needed
Metadata <- genetable_normdata %>%
  select(sample, type) %>%
  distinct()

ann_pcoa <- inner_join(pcoa_df, Metadata, by = "sample")

ann_pcoa_sorted <- ann_pcoa %>%
  distinct() %>%
  remove_rownames() %>%
  column_to_rownames(var = "sample")

nrow(clean_matrix)
nrow(ann_pcoa_sorted)

# PERMANOVA
permanova_results <- adonis2(
  bray_curtis_dist ~ type,
  data         = ann_pcoa_sorted,
  permutations = 9999,
  method       = "bray"
)
print(permanova_results)

r2_value        <- permanova_results$R2[1]
p_value         <- permanova_results$`Pr(>F)`[1]
permanova_label <- paste0("PERMANOVA R\u00B2: ", round(r2_value, 3),
                           "|", "P-value: ", round(p_value, 4))

# PERMDISP
betadisp <- betadisper(bray_curtis_dist, ann_pcoa_sorted$type)
permutest(betadisp)

# Colour palette
l_palette <- c("Exhale" = "#93e9f5", "Sputum" = "#fcf260")

# PCoA plot
pcaoa_plot <- ggplot(ann_pcoa, aes(x = PC1, y = PC2, label = sample, colour = type)) +
  geom_point(size = 3) +
  stat_ellipse(lwd = 0.8, linetype = "dashed") +
  scale_color_manual(values = l_palette) +
  labs(
    title = "",
    color = "Sample group",
    x     = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y     = paste0("PC2 (", round(var_explained[2], 2), "%)")
  ) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 12)) +
  annotate("text", x = Inf, y = Inf,
           label  = permanova_label,
           hjust  = 1.5, vjust = 0.9,
           size   = 4, color = "black")

pcaoa_plot
plotly::ggplotly(pcaoa_plot)

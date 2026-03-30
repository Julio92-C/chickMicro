# Purpose: PCoA of taxonomic profiles from the abricate+Bracken merged dataset,
#          using Bray-Curtis distances with Hellinger transformation, PERMANOVA,
#          and PERMDISP tests.
# Input:   data/abri_kraken2Bracken_merged.csv (cols: sample, name, sampleCount,
#                                               Treatment, DATABASE, ...)
# Output:  PCoA plot rendered to the active graphics device.

library(vegan)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(plotly)
library(paletteer)
library(funrar)

# Load dataset
abri_kraken2_merged <- read_csv("data/abri_kraken2Bracken_merged.csv")

# Build sample x taxon log-count matrix
taxaTreatment_dataset <- abri_kraken2_merged %>%
  filter(!grepl("Terrabacteria group|Bacteroidota/Chlorobiota group|FCB group", name)) %>%
  filter(str_count(name, "\\S+") > 1) %>%
  mutate(log_count = log(sampleCount + 1)) %>%
  select(sample, name, log_count) %>%
  group_by(sample) %>%
  tidyr::pivot_wider(
    names_from  = name,
    values_from = log_count,
    values_fill = 0,
    values_fn   = sum
  ) %>%
  distinct(sample, .keep_all = TRUE) %>%
  column_to_rownames("sample")

abundance_matrix <- as.matrix(taxaTreatment_dataset)

# Hellinger transformation and Bray-Curtis distance
relative_abundance_matrix <- decostand(abundance_matrix, method = "hellinger")

sum(is.na(relative_abundance_matrix))
sum(is.nan(as.matrix(relative_abundance_matrix)))

clean_matrix     <- na.omit(relative_abundance_matrix)
bray_curtis_dist <- vegdist(clean_matrix, method = "bray")

# PCoA
pcoa_result   <- cmdscale(bray_curtis_dist, eig = TRUE, k = 2)
eig_values    <- pcoa_result$eig
var_explained <- eig_values / sum(eig_values) * 100

pcoa_df <- data.frame(
  sample = rownames(clean_matrix),
  PC1    = pcoa_result$points[, 1],
  PC2    = pcoa_result$points[, 2]
)

Metadata <- abri_kraken2_merged %>%
  select(sample, Treatment) %>%
  distinct()

ann_pcoa <- inner_join(pcoa_df, Metadata, by = "sample")

ann_pcoa_sorted <- ann_pcoa %>%
  distinct() %>%
  remove_rownames() %>%
  column_to_rownames(var = "sample")

nrow(abundance_matrix)
nrow(ann_pcoa_sorted)

# PERMANOVA
permanova_results <- adonis2(
  bray_curtis_dist ~ Treatment,
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
betadisp <- betadisper(bray_curtis_dist, ann_pcoa_sorted$Treatment)
permutest(betadisp)

# Colour palette
l_palette <- c("Reference diet" = "#1f77b4",
               "Dulce"          = "#2ca02c",
               "Soyabean meal"  = "#ff7f0e")

# PCoA plot
pcaoa_plot <- ggplot(ann_pcoa, aes(x = PC1, y = PC2, label = sample, colour = Treatment)) +
  geom_point(size = 3) +
  stat_ellipse(lwd = 0.8, linetype = "dashed") +
  scale_color_manual(values = l_palette) +
  labs(
    title = "",
    color = "Treatment groups",
    x     = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y     = paste0("PC2 (", round(var_explained[2], 2), "%)")
  ) +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  annotate("text", x = Inf, y = Inf,
           label  = permanova_label,
           hjust  = 1.5, vjust = 0.9,
           size   = 4, color = "black")

pcaoa_plot
plotly::ggplotly(pcaoa_plot)

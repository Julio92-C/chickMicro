# Purpose: PCoA of ARG (CARD) gene profiles using Bray-Curtis distances, with
#          treatment group colouring and an optional plotly interactive output.
# Input:   data/abri_kraken2_merged.csv (cols: sample, name, GENE, Treatment, ...)
# Output:  PCoA plot rendered to the active graphics device; optionally save HTML
#          by uncommenting the saveWidget call below.

library(vegan)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(plotly)
library(paletteer)
library(htmlwidgets)

# Load dataset
abri_kraken2_merged <- read_csv("data/abri_kraken2_merged.csv")

# Prepare GENE x sample presence/absence matrix
abri_kraken2_arranged <- abri_kraken2_merged %>%
  arrange(GENE) %>%
  select(sample, name, GENE, Treatment)

unique_gene <- length(unique(abri_kraken2_arranged$GENE))

df_wide_ann_rows <- abri_kraken2_arranged %>%
  arrange(sample) %>%
  pivot_wider(names_from = sample, values_from = c(name), values_fn = length) %>%
  group_by(GENE) %>%
  summarise(across(2:19, ~ paste(., collapse = ", "))) %>%
  mutate(across(2:19, ~ str_replace_all(., "(NA,|,NA|,NA,|NA|,| )", "")))

df_wide_ann_rows <- df_wide_ann_rows %>%
  na.omit() %>%
  remove_rownames() %>%
  column_to_rownames(var = "GENE") %>%
  mutate(across(everything(), ~ if_else(. == "", NA_character_, .)))

df_clean <- df_wide_ann_rows[
  !apply(df_wide_ann_rows, 1, function(row) all(row == "" | is.na(row))), ]

sample2arg_data_matrix <- data.matrix(df_clean, rownames.force = NA)
sum(is.na(sample2arg_data_matrix))
sample2arg_data_matrix <- sample2arg_data_matrix %>% replace(is.na(.), 0)

# Bray-Curtis PCoA
bray_curtis_dist <- vegdist(sample2arg_data_matrix, method = "bray")
pcoa_result      <- cmdscale(bray_curtis_dist, eig = TRUE, k = 2)
eig_values       <- pcoa_result$eig
var_explained    <- eig_values / sum(eig_values) * 100

pcoa_df <- data.frame(
  GENE = rownames(sample2arg_data_matrix),
  PC1  = pcoa_result$points[, 1],
  PC2  = pcoa_result$points[, 2]
)

treatmentData <- abri_kraken2_arranged %>%
  select(sample, GENE, Treatment) %>%
  distinct()

ann_pcoa <- inner_join(pcoa_df, treatmentData, by = "GENE")

treatment_col <- as.factor(ann_pcoa$Treatment)
l             <- length(levels(treatment_col))
l_palette     <- paletteer_d("ggsci::default_igv", n = l)

# PCoA plot
pcaoa_plot <- ggplot(ann_pcoa, aes(x = PC1, y = PC2, colour = Treatment)) +
  geom_point(size = 3) +
  stat_ellipse(type = "norm", lwd = 0.8) +
  scale_color_manual(values = l_palette) +
  labs(
    title  = "",
    color  = "Treatment",
    x      = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y      = paste0("PC2 (", round(var_explained[2], 2), "%)")
  ) +
  theme_classic() +
  theme(text = element_text(size = 16)) +
  scale_x_continuous(limits = c(-0.8, NA)) +
  scale_y_continuous(limits = c(-0.6, NA))

pcaoa_plot

fig <- plotly::ggplotly(pcaoa_plot)
fig

# Uncomment to save interactive plot and data:
# saveWidget(widget = fig, file = "output/PCoA_plot.html", selfcontained = TRUE)
# write.csv(ann_pcoa, "output/ann_pcoa.csv", row.names = FALSE)

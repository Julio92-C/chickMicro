# Purpose: Binary presence/absence pheatmap of taxa carrying AMR, VFs, and MGEs,
#          annotated by treatment group, drug class, resistance type, and MGE burden.
# Input:   data/abri_kraken2Bracken_merged.csv (cols: sample, name, GENE, PRODUCT,
#                                               RESISTANCE, DATABASE, Treatment, sampleCount, ...)
# Output:  Heatmap rendered to the active graphics device.

library(readr)
library(dplyr)
library(tidyverse)
library(paletteer)
library(reshape2)
library(pheatmap)
library(stringr)

# Load dataset
abri_kraken2_clean <- read_csv("data/abri_kraken2Bracken_merged.csv")

# Filter to species-level taxa with sufficient sample coverage
abri_kraken2_clean <- abri_kraken2_clean %>%
  filter(!grepl("Terrabacteria group|Bacteroidota/Chlorobiota group|FCB group", name)) %>%
  filter(!(str_count(name, "\\S+") == 1 & name != "Enterococcus")) %>%
  mutate(
    name = str_replace(name, "^(\\w)\\w+\\s(\\w+)", "\\1. \\2"),
    name = if_else(
      str_length(name) > 30,
      str_replace(name, "^((?:\\S+\\s+){2}).*", "\\1"),
      name
    )
  ) %>%
  mutate(name = case_when(
    grepl("\\[Clostridium\\] asparagiforme ", name) ~ "C. asparagiforme",
    grepl("\\[Clostridium\\] scindens", name)       ~ "C. scindens",
    grepl("\\[Ruminococcus\\] lactaris", name)      ~ "R. lactaris ",
    grepl("\\[Ruminococcus\\] torques", name)       ~ "R. torques",
    grepl("\\[Clostridium] innocuum", name)         ~ "C. innocuum",
    TRUE ~ name
  ),
  name = ifelse(grepl("u. Subdoligranulum", name), "U. Subdoligranulum", name)) %>%
  filter(sampleCount > 10)

uniqueTaxa <- length(unique(abri_kraken2_clean$name))

# Wide binary matrix: taxa x samples
df_wide <- abri_kraken2_clean %>%
  arrange(Treatment) %>%
  select(sample, taxid, name, GENE, Treatment) %>%
  pivot_wider(names_from = sample, values_from = c(taxid), values_fn = length) %>%
  group_by(name) %>%
  summarise(across(3:20, ~ paste(., collapse = ", "))) %>%
  mutate(across(2:19, ~ str_replace_all(., "(NA,|,NA|,NA,|NA|,| )", "")))

processed_mgs2arg_data <- df_wide %>%
  remove_rownames() %>%
  column_to_rownames(var = "name") %>%
  mutate(across(everything(), ~ na_if(., "")))

sample2arg_data <- processed_mgs2arg_data %>%
  filter(rowSums(is.na(.)) < ncol(.))

sample2arg_data_matrix <- data.matrix(sample2arg_data, rownames.force = NA) %>%
  replace(is.na(.), 0)

sum(is.na(sample2arg_data_matrix))

# Convert to binary presence/absence
sample2arg_data_Bmatrix <- as.matrix((sample2arg_data_matrix > 0) + 0)

# Sort rows by decreasing row sum
rearrange_matrix <- function(mat) {
  mat[order(rowSums(mat), decreasing = TRUE), ]
}
sorted_matrix <- rearrange_matrix(sample2arg_data_Bmatrix)

# Clustering callback for SVD-guided dendrogram ordering
callback <- function(hc, mat) {
  sv   <- svd(t(mat))$v[, 1]
  dend <- reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

# Column annotation: treatment group
sample2arg_data_transposed <- as.data.frame(t(sample2arg_data)) %>%
  mutate(sample = rownames(.)) %>%
  select(sample)

MetadataLocations <- abri_kraken2_clean %>%
  select(sample, Treatment) %>%
  distinct()

ann_col <- inner_join(sample2arg_data_transposed, MetadataLocations, by = "sample") %>%
  na.omit() %>%
  remove_rownames() %>%
  column_to_rownames(var = "sample")

d_palette  <- paletteer_d("ggsci::category10_d3", n = nlevels(as.factor(ann_col$Treatment)))
df4        <- ann_col %>% distinct(Treatment) %>% mutate(color = d_palette)
type_colors <- setNames(as.character(df4$color), df4$Treatment)

# Row annotation: VFs / AMR / MGE classification
df_wide_ann_rows <- abri_kraken2_clean %>%
  arrange(sample) %>%
  filter(grepl("card|vfdb|plasmidfinder", DATABASE)) %>%
  pivot_wider(names_from = DATABASE, values_from = GENE, values_fn = list) %>%
  group_by(name) %>%
  summarise(across(21:23, ~ paste(., collapse = ", "))) %>%
  mutate(across(2:4, ~ str_replace_all(., "(NULL,|,NULL|,NULL,|NULL| )", ""))) %>%
  na.omit() %>%
  remove_rownames() %>%
  column_to_rownames(var = "name") %>%
  mutate(across(everything(), ~ if_else(. == "", NA_character_, .)))

classify <- function(input_string, type) {
  if (is.null(input_string) || is.na(input_string)) return(NA)
  elements  <- unlist(strsplit(input_string, ","))
  elements  <- gsub("^\\s+|\\s+$", "", elements)
  elements  <- elements[elements != ""]
  n_unique  <- length(unique(elements))
  if (type == "virulence") return(ifelse(n_unique == 1, "Virulent", "Hypervirulent"))
  if (type == "amr")       return(ifelse(n_unique == 1, "Drug-resistant", "Multidrug-resistant"))
  if (type == "MGE")       return(ifelse(n_unique == 1, "Single MGE", "Multi-MGEs"))
  return(NA)
}

ann_rows <- df_wide_ann_rows %>%
  mutate(VFs  = sapply(vfdb,         classify, type = "virulence"),
         AMR  = sapply(card,          classify, type = "amr"),
         MGEs = sapply(plasmidfinder, classify, type = "MGE")) %>%
  select(VFs, AMR, MGEs)

# Drug resistance annotation
clean_string <- function(x) str_to_title(str_replace_all(x, "_", " "))

Resistance_category <- abri_kraken2_clean %>%
  distinct(name, GENE, RESISTANCE) %>%
  mutate(RESISTANCE = sapply(RESISTANCE, clean_string)) %>%
  group_by(name) %>%
  summarise(RESISTANCE = paste(unique(RESISTANCE), collapse = ", ")) %>%
  column_to_rownames(var = "name") %>%
  mutate(across(1, ~ str_replace_all(., "(NA,|,NA|,NA,|NA|,| )", ""))) %>%
  mutate(RESISTANCE = ifelse(nchar(RESISTANCE) > 20, "Multi-drug", RESISTANCE)) %>%
  mutate(across(everything(), ~ if_else(. == "", NA_character_, .))) %>%
  select(RESISTANCE)

ann_row_all <- merge(ann_rows, Resistance_category, by = 0) %>%
  rename(DRUG = RESISTANCE) %>%
  column_to_rownames(var = "Row.names") %>%
  select(DRUG, AMR, VFs, MGEs)

r_palette <- paletteer_d("ggsci::default_igv", n = length(unique(ann_row_all$DRUG)))

ann_row_drug_col <- ann_row_all %>%
  select(DRUG) %>%
  na.omit() %>%
  distinct(DRUG) %>%
  mutate(color = r_palette)

drug_colors <- setNames(as.character(ann_row_drug_col$color), ann_row_drug_col$DRUG)

ann_colors <- list(
  Treatment = c("Reference diet" = "#1f77b4",
                "Dulce"          = "#2ca02c",
                "Soyabean meal"  = "#ff7f0e"),
  DRUG = drug_colors,
  AMR  = c(`Drug-resistant` = "#a02d02", `Multidrug-resistant` = "#522501"),
  MGEs = c(`Multi-MGEs` = "#5c0404", `Single MGE` = "#fcd2d2"),
  VFs  = c(`Virulent` = "#fce5d2", `Hypervirulent` = "#c47a02")
)

# Pheatmap
heatmap_plot <- pheatmap(
  sorted_matrix,
  display_numbers      = FALSE,
  cluster_cols         = FALSE,
  cluster_rows         = FALSE,
  scale                = "none",
  clustering_callback  = callback,
  border_color         = "NA",
  color                = c("#CCCCCCFF", "#666666FF"),
  legend_breaks        = c(0, 1),
  legend_labels        = c("Absent", "Present"),
  annotation_row       = ann_row_all,
  annotation_col       = ann_col[1],
  show_rownames        = TRUE,
  annotation_colors    = ann_colors,
  fontsize_row         = 14,
  fontsize_col         = 14,
  fontsize             = 14
)

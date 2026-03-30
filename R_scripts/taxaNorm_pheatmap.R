# Load the libraries
library(readr)
library(dplyr)
library(tidyverse)
library(paletteer)
library(reshape2)
library(pheatmap)
library(stringr)
library(funrar)
library(vegan)


# Import Bracken arrange report 
bracken_arranged <- read_csv("../data/bracken_arranged.csv")
#View(bracken_arranged)
colnames(bracken_arranged)

# Load metadata
chicken_metadata <- read_csv("../data/chicken_metadata.csv")
chicken_metadata <- na.omit(chicken_metadata)
View(chicken_metadata)  

# Load candidates from DA analysis
aldex_top_candidates_summary <- read_csv("../data/aldex_top_candidates_summary.csv")
View(aldex_top_candidates_summary)

# Modify strings function
modify_strings <- function(df, column_name, target_string, replacement_string) {
  df <- df %>%
    mutate(!!sym(column_name) := str_replace_all(!!sym(column_name), fixed(target_string), replacement_string))
  return(df)
}


# Call Modify String function
chicken_metadata <- modify_strings(chicken_metadata, "Treatment", 
                                   "Reference diet (Wheat -soyabean)", 
                                   "Reference diet")
chicken_metadata <- modify_strings(chicken_metadata, "Treatment", 
                                   "Seaweed", 
                                   "Dulce")

# Pivot longer Bracken report
bracken_pivoted <- bracken_arranged %>%
  pivot_longer(cols = 4:21, names_to = "sample", values_to = "count") %>%
  filter(!grepl("root|Homo sapiens|cellular organisms|unclassified|Bacteria|environmental samples", name)) %>%
  filter(!is.na(name)) %>%
  group_by(sample) %>%
  #filter(count > 5) %>%
  mutate(Percentage = (count / sum(count))*100) %>%
  ungroup()

# Merge bracken_arranged and metadata dataframe
bracken_merged <- merge(bracken_pivoted, chicken_metadata, by ="sample")

# Candidates list
candidates_list <- aldex_top_candidates_summary$taxon

# Wheat-soyabean dataset
taxaTreatment_dataset <- bracken_merged %>%
  filter(!grepl("Terrabacteria group|Bacteroidota/Chlorobiota group|FCB group|Cyanobacteriota/Melainabacteria group|Rhizobium/Agrobacterium group|Candida/Lodderomyces clade|Sinorhizobium/Ensifer group|Plasmodium", name)) %>%
  filter(str_count(name, "\\S+") > 1) %>%
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
    grepl("\\[Clostridium\\] scindens", name) ~ "C. scindens",
    grepl("\\[Ruminococcus\\] lactaris", name) ~ "R. lactaris ",
    grepl("\\[Ruminococcus\\] torques", name) ~ "R. torques",
    grepl("\\[Clostridium\\] innocuum", name) ~ "C. innocuum",
    grepl("\\[Clostridium\\] colinum", name) ~ "C. colinum",
    grepl("\\[Clostridium\\] hylemonae", name) ~ "C. hylemonae",
    TRUE ~ name
  ),
  name = ifelse(grepl("u. Subdoligranulum", name), "U. Subdoligranulum", name)) %>%
  arrange(Treatment) %>%  
  select(sample, name, Percentage) %>%
  group_by(sample) %>%
    tidyr::pivot_wider(
    names_from = sample,
    values_from = Percentage,
    values_fill = 0,
    values_fn = sum # Use `sum` to aggregate duplicate values
  ) %>%
  distinct(name, .keep_all = TRUE) %>%
  filter(name %in% candidates_list) %>%
  column_to_rownames("name")

# Convert the data frame to a matrix
abundance_matrix <- as.matrix(taxaTreatment_dataset)

# Hellinger transformation
#relative_abundance_matrix <- decostand(abundance_matrix, method = "hellinger")


# Check that scaling worked
range(abundance_matrix)     # Should show 0 to 1
rowSums(abundance_matrix)   # Should all be 1

# Rearranging a Matrix by Row Count
rearrange_matrix <- function(mat) {
  # Calculate the row sums (count of '1's)
  ordered_indices <- order(rowSums(mat), decreasing = TRUE)
  
  # Rearrange the matrix
  return(mat[ordered_indices, ])
}

# Rearranging the gene matrix
taxaTable_matrix_sorted <- rearrange_matrix(abundance_matrix)

# Modify ordering of the clusters using clustering callback option
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}


# Annotations col names
# Transpose sample2arg_data
taxaTable_filtered_transposed <- as.data.frame(t(taxaTreatment_dataset))

# Convert rownames into column and select sample
taxaTable_filtered_transposed <- taxaTable_filtered_transposed %>%
  mutate(sample = rownames(.)) %>%
  select(sample)

# Select sample and treatmet category
Metadata_Treatment <- bracken_merged  %>%
  select(sample,Treatment) %>%
  distinct()

# Merge the metadata and sample dataset
ann_col <- inner_join(taxaTable_filtered_transposed, Metadata_Treatment, by = "sample") %>%
  na.omit() %>%
  remove_rownames() %>%
  column_to_rownames(var = "sample")

# Treatment dataset
taxaTreatment_arraged <- bracken_merged %>%
  filter(!grepl("Terrabacteria group|Bacteroidota/Chlorobiota group|FCB group|Cyanobacteriota/Melainabacteria group|Rhizobium/Agrobacterium group|Candida/Lodderomyces clade|Sinorhizobium/Ensifer group|Plasmodium", name)) %>%
  filter(str_count(name, "\\S+") > 1) %>%
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
    grepl("\\[Clostridium\\] scindens", name) ~ "C. scindens",
    grepl("\\[Ruminococcus\\] lactaris", name) ~ "R. lactaris ",
    grepl("\\[Ruminococcus\\] torques", name) ~ "R. torques",
    grepl("\\[Clostridium\\] innocuum", name) ~ "C. innocuum",
    grepl("\\[Clostridium\\] colinum", name) ~ "C. colinum",
    grepl("\\[Clostridium\\] hylemonae", name) ~ "C. hylemonae",
    TRUE ~ name
  ),
  name = ifelse(grepl("u. Subdoligranulum", name), "U. Subdoligranulum", name)) %>%
  select(name, Treatment, log_count) %>%
  tidyr::pivot_wider(
    names_from = Treatment,
    values_from = log_count,
    values_fill = 0,
    values_fn = sum # Use `sum` to aggregate duplicate values
  ) %>%
  distinct(name, .keep_all = TRUE) %>%
  column_to_rownames("name")

## Convert abundance to binary presence/absence
presence_matrix <- as.matrix((taxaTreatment_arraged > 0))

# Count how many treatments each taxon appears in
presence_count <- rowSums(presence_matrix)

# Create dataframe with taxa + category
taxa_category <- data.frame(
  Taxon = rownames(presence_matrix),
  Category = case_when(
    presence_count == 3 ~ "Core", # High occupancy: consistent across all samples
    presence_count == 2 ~ "Accessory", # Moderate occupancy: intermediate distribution
    presence_count == 1 ~ "Unique", # Low occupancy: sample-specific or rare
    TRUE ~ "Absent"
  )
)

# Data categories
#write.csv(taxa_category, "../data/taxa_category.csv", row.names = F)

# Set the row names
taxa_category <- taxa_category %>%
  remove_rownames() %>%
  column_to_rownames(var = "Taxon")

# Print counts for each category
core_n       <- sum(taxa_category$Category == "Core")
accessory_n  <- sum(taxa_category$Category == "Accessory")
unique_n     <- sum(taxa_category$Category == "Unique")

# Compute gaps automatically
gaps_row <- c(core_n, core_n + accessory_n)

# Compute gaps automatically
gaps_col <- c(6, 12)

# Annotation colors customization
ann_colors <- list(
  Treatment = c("Reference diet" = "#1f77b4",
                "Dulce" = "#2ca02c",
                "Soyabean meal" = "#ff7f0e"),
  Category = c("Core" = "#15b379", 
               "Accessory" = "#b9e9ed",
               "Unique" = "#f2615a")
  
)


#Create heatmap using pheatmap package ##
heatmap_plot <- pheatmap(taxaTable_matrix_sorted, display_numbers = FALSE, cluster_cols = F, cluster_rows = F,
                         scale = "none",
                         #clustering_callback = callback,
                         color = colorRampPalette(c("#15b379","yellow", "#f2615a"))(10),
                         #legend_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                         #legend_labels = c("0", "0.25", "0.5", "0.75", "1"),
                         border_color = NA,
                         show_rownames = T,
                         # cutree_cols = 15,
                         # cutree_rows = 20,
                         #annotation_row = taxa_category,
                         #gaps_row = gaps_row,
                         gaps_col = gaps_col,
                         annotation_col = ann_col,
                         annotation_colors = ann_colors,
                         fontsize_row = 14,  # Adjust this value for row names
                         fontsize_col = 14,   # Adjust this value for column names
                         fontsize = 14  # Adjust this value for the legend text
                         
                         
)

# Plot pheatmap
heatmap_plot


# Clear environment
rm(list = ls())

# Clear packages
detach("package:datasets", unload = TRUE)

# Clear plots
dev.off()  # But only if there IS a plot

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)


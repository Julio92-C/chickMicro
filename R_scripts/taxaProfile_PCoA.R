
# Load libraries
library(vegan)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(plotly)
library(paletteer)
library(htmlwidgets)


# Import Bracken arrange report 
bracken_arranged <- read_csv("../data/bracken_arranged.csv")
View(bracken_arranged)
colnames(bracken_arranged)

# Load metadata
chicken_metadata <- read_csv("../data/chicken_metadata.csv")
chicken_metadata <- na.omit(chicken_metadata)
View(chicken_metadata)  


# Pivot longer Bracken report
bracken_pivoted <- bracken_arranged %>%
  pivot_longer(cols = 4:21, names_to = "sample", values_to = "count") %>%
  filter(!grepl("root|Homo sapiens|cellular organisms|unclassified|Bacteria|environmental samples", name)) %>%
  filter(!is.na(name)) %>%
  group_by(sample) %>%
  filter(count > 5) %>%
  mutate(Percentage = (count / sum(count)) * 100) %>%
  ungroup()

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

# Merge bracken_arranged and metadata dataframe
bracken_merged <- merge(bracken_pivoted, chicken_metadata, by ="sample")

# Treatment categories
treat_categories <- unique(bracken_merged$Treatment)


# Arrange dataset by gene variable
abri_kraken2_arranged <- bracken_merged %>%
  arrange(name) %>%
  select(sample, name, taxid, Treatment)

# Count unique number of taxa
unique_taxa <- length(unique(abri_kraken2_arranged$name))

# Annotation rows
# Pivot the taxa values into columns
df_wide_ann_rows <- abri_kraken2_arranged %>%
  arrange(sample) %>%
  pivot_wider(names_from = sample, values_from = c(taxid), values_fn = length) %>%
  group_by(name) %>%
  summarise(across(everything(), ~ paste(., collapse = ", "))) %>%
  mutate(across(-name, ~ str_replace_all(., "(NA,|,NA|,NA,|NA|,| )", "")))

# Assign the values of GENE as row names and empty rows as NA value
df_wide_ann_rows <- df_wide_ann_rows %>%
  na.omit()  %>%
  remove_rownames %>%
  column_to_rownames(var="name") %>%
  mutate(across(everything(), ~ if_else(. == "", NA_character_, .)))


# Remove rows where all values are NA or empty strings
df_clean <- df_wide_ann_rows[!apply(df_wide_ann_rows, 1, function(row) all(row == "" | is.na(row))), ]


# Convert from a dataframe to a matrix
sample2arg_data_matrix <- data.matrix(df_clean, rownames.force = NA)

# Count NA in data sets
sum(is.na(sample2arg_data_matrix))

# Replace all NA values with 0 values
sample2arg_data_matrix <- sample2arg_data_matrix %>% replace(is.na(.), 0)

# Calculate Bray-Curtis distances
bray_curtis_dist <- vegdist(sample2arg_data_matrix, method = "bray")

# Perform PCoA
pcoa_result <- cmdscale(bray_curtis_dist, eig = TRUE, k = 2)

# Calculate percentage of variance explained
eig_values <- pcoa_result$eig
var_explained <- eig_values / sum(eig_values) * 100

# Create a data frame for plotting
pcoa_df <- data.frame(name = rownames(sample2arg_data_matrix),
                      PC1 = pcoa_result$points[, 1],
                      PC2 = pcoa_result$points[, 2])


# Metadata dataframe
treatmentData <- abri_kraken2_arranged %>%
  select(sample, name, Treatment) %>%
  distinct()

# Merge the metadata and sample dataset
# ann_pcoa <- inner_join(pcoa_df, MetadataLocations, by="sample")
ann_pcoa <- inner_join(pcoa_df, treatmentData, by="name")

# Annotation col colors for Treatment
# Description variable
treatment_col <- as.factor(ann_pcoa$Treatment)

# Generate a color palette based on the number of levels in gene_family
l <- length(levels(treatment_col))
l_palette <- paletteer_d("ggsci::default_igv", n = l)
l_palette

# Plot the results
pcaoa_plot <- ggplot(ann_pcoa, aes(x = PC1, y = PC2, colour = Treatment)) +
  geom_point(size = 3) +
  #geom_text(aes(label = GENE), vjust = 1.5, hjust = 1.5) +
  # stat_ellipse(geom = "polygon", aes(fill = Treatment), alpha = 0.2) +
  stat_ellipse(type = "norm", lwd = 0.8) + # Change line width
  # geom_jitter(width = 0.02, height = 0.02) +
  scale_color_manual(values = l_palette) +
  #scale_fill_manual(values = l_palette) +
  labs(title = "",
       color = "Treatment",
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 2), "%)")
  ) +
  theme_classic() +
  theme(text = element_text(size = 16)) +
  scale_x_continuous(limits = c(-0.6, NA)) +  # Set x-axis to start at -0.5
  scale_y_continuous(limits = c(-0.6, NA))    # Set y-axis to start at -0.5

# ggplot graph
pcaoa_plot

# ploty graph
fig <- plotly::ggplotly(pcaoa_plot)
fig

# Same to a html plot
saveWidget(widget = fig, file = "../figures/PCoA_plot.html", selfcontained = TRUE)

#   Save ann_pcoa to a csv file
write.csv(ann_pcoa, "../data/ann_pcoa.csv", row.names = FALSE)

# CLEAN UP #################################################

# Clear environment
rm(list = ls())

# Clear packages
detach("package:datasets", unload = TRUE)

# Clear plots
dev.off()  # But only if there IS a plot

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)
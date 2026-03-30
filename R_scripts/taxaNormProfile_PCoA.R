# Load libraries
library(vegan)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(plotly)
library(paletteer)
library(funrar)


# Import Bracken arrange report 
bracken_arranged <- read_csv("../data/bracken_arranged.csv")
#View(bracken_arranged)
colnames(bracken_arranged)

# Load metadata
chicken_metadata <- read_csv("../data/chicken_metadata.csv")
chicken_metadata <- na.omit(chicken_metadata)
View(chicken_metadata)  

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
  filter(count > 5) %>%
  mutate(log_count = log(count + 1)) %>%
  ungroup()

# Merge bracken_arranged and metadata dataframe
bracken_merged <- merge(bracken_pivoted, chicken_metadata, by ="sample")


# Wheat-soyabean dataset
taxaTreatment_dataset <- bracken_merged %>%
  filter(!grepl("Terrabacteria group|Bacteroidota/Chlorobiota group|FCB group", name)) %>%
  filter(str_count(name, "\\S+") > 1) %>%
  select(sample, name, log_count) %>%
  group_by(sample) %>%
  tidyr::pivot_wider(
    names_from = name,
    values_from = log_count,
    values_fill = 0,
    values_fn = sum # Use `sum` to aggregate duplicate values
  ) %>%
  distinct(sample, .keep_all = TRUE) %>%
  column_to_rownames("sample")
  
  
# Convert the data frame to a matrix
abundance_matrix <- as.matrix(taxaTreatment_dataset)

# Hellinger transformation
relative_abundance_matrix <- decostand(abundance_matrix, method = "hellinger")


# Check how many missing values you have
sum(is.na(relative_abundance_matrix))
sum(is.nan(as.matrix(relative_abundance_matrix))) # Check for NaNs as well

# Remove rows with any missing values
clean_matrix <- na.omit(relative_abundance_matrix)

# Compute Bray-Curtis distance
bray_curtis_dist <- vegdist(clean_matrix, method = "bray")

# Perform PCoA
pcoa_result <- cmdscale(bray_curtis_dist, eig = TRUE, k = 2)

# Calculate percentage of variance explained
eig_values <- pcoa_result$eig
var_explained <- eig_values / sum(eig_values) * 100

# Create a data frame for plotting
pcoa_df <- data.frame(sample = rownames(clean_matrix),
                      PC1 = pcoa_result$points[, 1],
                      PC2 = pcoa_result$points[, 2])


# Metadata subset
Metadata <- bracken_merged %>%
  select(sample, Treatment) %>%
  distinct()

# # Merge the metadata and sample dataset
ann_pcoa <- inner_join(pcoa_df, Metadata, by="sample")

# Re-arrange ann_pcoa
ann_pcoa_sorted <- ann_pcoa %>%
  distinct() %>%
  remove_rownames %>%
  column_to_rownames(var="sample")

# Check that the number of row match
nrow(abundance_matrix) # Should be the number of samples you have
nrow(ann_pcoa_sorted)     # Should be the number of samples you have


# Run the PERMANOVA test
permanova_results <- adonis2(
  bray_curtis_dist ~ Treatment,          # Formula: distance matrix explained by 'Group'
  data = ann_pcoa_sorted,           # Data frame where 'Group' variable is found
  permutations = 9999,           # Number of permutations (999 or 9999 for publication)
  method = "bray"               # Specify distance method if not already specified in the matrix
)

# Print the results
print(permanova_results)


# Extract R2 and P-value for the 'type' variable
# The 'type' variable is on the first row of results (after the 'Df' column, before 'Residual')
# The results are stored within the data frame structure of the permanova_results object
r2_value <- permanova_results$R2[1]
p_value <- permanova_results$`Pr(>F)`[1]

# Format the results into a single label string for the plot
# Use paste0 to combine text and round the values
permanova_label <- paste0(
  "PERMANOVA R\u00B2: ", round(r2_value, 3), "|", "P-value: ", round(p_value, 4)
)

# Test for homogeneity of dispersion (PERMDISP)
# PERMANOVA assumes equal dispersion across groups.
betadisp <- betadisper(bray_curtis_dist, ann_pcoa_sorted$Treatment)
permutest(betadisp)


# Annotation col colors for Location
l_palette <- c("Reference diet" = "#1f77b4",
               "Dulce" = "#2ca02c",
               "Soyabean meal" = "#ff7f0e")

# Plot the results
pcaoa_plot <- ggplot(ann_pcoa, aes(x = PC1, y = PC2, label = sample, colour = Treatment)) +
  geom_point(size = 3) +
  #geom_text(vjust = 1.5, hjust = 1.5) +
  #stat_ellipse(geom = "polygon", aes(fill = location), alpha = 0.2) +
  stat_ellipse(lwd = 0.8, # Change line width
               linetype = "dashed") + # Use a dashed line style
  # Other options include "dotted", "dotdash", "longdash", "twodash" 
  scale_color_manual(values = l_palette) +
  #scale_fill_manual(values = l_palette) +
  #geom_path(aes(group = Treatment_Groups)) +
  #coord_cartesian(xlim = c(-0.6, 1), ylim = c(-0.5, 0.3)) + # Zoom the plot area
  labs(title = "",
       color = "Treatment groups",
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 2), "%)")
  ) +
  theme_classic() +
  theme(text = element_text(size = 14)) +
  # --- 3. Add the PERMANOVA results as text annotation ---
  annotate(
    "text", 
    x = Inf, y = Inf,             # Position the text at the top right corner
    label = permanova_label,      # The string we created earlier
    hjust = 1.5, vjust = 0.9,     # Adjust justification slightly inwards
    size = 4,                     # Control font size
    color = "black"               # Control font color
  )

# ggplot graph
pcaoa_plot

plotly::ggplotly(pcaoa_plot)



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



# Load libraries
library(ggplot2)
library(gggenes)

# Load dataset
abri_kraken2_merged <- read_csv("../data/abri_kraken2_merged.csv")


# Filter the dataset by Escherichia species and VFs genes
df_species <- abri_kraken2_merged %>%
  filter(grepl("Escherichia", name)) %>%
  #filter(grepl("Enterobacteriaceae", name)) %>%
  #filter(grepl("Clostridioides difficile", name)) %>%
  filter(grepl("vfdb", DATABASE)) %>%
  #filter(DATABASE == "card") %>%
  #filter(name != "Clostridia") %>%
  #filter(GENE != "AAC(6')-Ie-APH(2'')-Ia_bifunctional_protein") %>%
  # filter(DATABASE != "ncbi") %>%
  distinct()

# Unique genes
uniqueGene <- length(unique(df_species$GENE))
uniqueSpecies <- length(unique(df_species$name))
colnames(df_species)

# Combine sample and name into one column
df_species_combined <- df_species %>%
  unite("genome", c(18,1), sep = " ") %>%
  group_by(Treatment) %>%
  distinct(genome, GENE, .keep_all = TRUE) 


# Remove NA values
# df_species_combined <- na.omit(df_species_combined)


# dummies <- make_alignment_dummies(
#   df_species_combined,
#   ggplot2::aes(xmin = START, xmax = END, y = genome, id = GENE),
#   on = "APH(3')-IIIa"
# )


# Annotation col colors for taxa
# Description variable
geneUnique <- as.factor(df_species_combined$GENE)

# Generate a color palette based on the number of levels in gene_family
t <- length(levels(geneUnique))
t_palette <- paletteer_d("ggsci::default_igv", n = t)
t_palette

#df <- palettes_d_names

# # Create the plot
p <- ggplot(df_species_combined, aes(xmin = START, xmax = END, y = genome, fill = GENE)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  #geom_blank(data = dummies) +
  #facet_wrap(~ genome, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3") +
  scale_fill_manual(values = t_palette) +
  theme_genes() +
  #geom_gene_label(aes(label = GENE), align = "centre", size = 6) +
  facet_grid(Treatment ~ ., scales = "free_y", switch = "y")

# Display the plot
print(p)




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

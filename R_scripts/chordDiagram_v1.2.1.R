# Load the libraries
library(readr)
library(circlize)
library(plotly)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(grid)
library(paletteer)

# Load dataset
abri_kraken2_clean <- read_csv("../data/abri_kraken2_filtered.csv")

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

# Merge bracken_arranged and metadata dataframe
abri_kraken2_merged <- merge(abri_kraken2_clean, chicken_metadata, by ="sample")

# Define modify_strings function
modify_strings <- function(df, column_name, target_strings, replacement_strings) {
  if (length(target_strings) != length(replacement_strings)) {
    stop("target_strings and replacement_strings must be of the same length")
  }
  for (i in seq_along(target_strings)) {
    df <- df %>%
      mutate(!!sym(column_name) := str_replace_all(!!sym(column_name), target_strings[i], replacement_strings[i]))
  }
  return(df)
}

# List of target and replacement strings
target_strings <- c("carbapenem;cephalosporin;penam")
replacement_strings <- c("Beta-lactam")

# Call Modify String function
abri_kraken2_filtered <- modify_strings(abri_kraken2_merged, "RESISTANCE", target_strings, replacement_strings)





# Filter and mutate the dataset
abri_kraken2_filtered <- abri_kraken2_filtered %>%
  filter(grepl("card", DATABASE)) %>%
  mutate(GENE = case_when(
    grepl("vanW_gene_in_vanB_cluster", GENE) ~ "vanWB",
    grepl("vanR_gene_in_vanB_cluster", GENE) ~ "vanRB",
    grepl("vanX_gene_in_vanB_cluster", GENE) ~ "vanXB",
    grepl("vanY_gene_in_vanB_cluster", GENE) ~ "vanYB",
    grepl("vanS_gene_in_vanB_cluster", GENE) ~ "vanSB",
    grepl("vanH_gene_in_vanB_cluster", GENE) ~ "vanHB",
    grepl("Escherichia_coli_ampC_beta-lactamase", GENE) ~ "ampC",
    grepl("Escherichia_coli_mdfA", GENE) ~ "mdfA",
    grepl("Escherichia_coli_emrE", GENE) ~ "emrE",
    TRUE ~ GENE
  ),
  GENE = ifelse(grepl("AAC(6')-Ie-APH(2'')-Ia_bifunctional_protein", GENE), "aac(6')-Ie/aph(2'')-Ia", GENE)) %>%
  mutate(name = ifelse(grepl("Chryseobacterium panacisoli", name), "C. panacisoli", name)) %>%
  mutate(RESISTANCE = ifelse(nchar(RESISTANCE) > 20, "Multi-drug", RESISTANCE)) %>%
  arrange(name)

# Clean the RESISTANCE strings
abri_kraken2_filtered <- abri_kraken2_filtered %>%
  mutate(RESISTANCE = str_replace_all(RESISTANCE, "_", " ") %>% str_to_title())

# Resistance gene category
pathogen_AMR <- abri_kraken2_filtered %>%
  filter(!grepl("Alistipes|Azorhizobium|cellular organisms|Bacillota|Clostridia|Lachnospiraceae|Bacteria|root|Enterobacterales|Enterobacteriaceae|Terrabacteria group|Homo sapiens|Pseudomonadota|Bacteroidales|Bifidobacterium|Bacteroidota", name)) %>%
  filter(!(str_count(name, "\\S+") == 1 & name != "Enterococcus")) %>%
  filter(Treatment == "Reference diet") %>%
  # filter(Treatment == "Soyabean meal") %>%
  #filter(Treatment == "Seaweed") %>%
  select(sample, name, GENE, RESISTANCE) %>%
  mutate(name = str_replace(name, "^(\\w)\\w+\\s(\\w+)", "\\1. \\2")) %>%
  distinct()

# Categories variables
name_category <- unique(as.factor(pathogen_AMR$name))
sample_category <- unique(as.factor(pathogen_AMR$sample))
gene_category <- unique(as.factor(pathogen_AMR$GENE))
resistant_category <- unique(as.factor(pathogen_AMR$RESISTANCE))

#pathogen_AMR[26, 2] = "aac(6')-Ie/aph(2'')-Ia"
# Reference diet
pathogen_AMR[64, 2] = "C. innocuum"
pathogen_AMR[65, 2] = "C. scindens"
pathogen_AMR[67:68, 2] = "u. Subdoligranulum"
# Soyabean meal
pathogen_AMR[59, 2] = "C. scindens"
pathogen_AMR[60, 2] = "C. scindens"
pathogen_AMR[61, 2] = "R torques"
pathogen_AMR[10, 3] = "aac(6')-Ie/aph(2'')-Ia"
pathogen_AMR[47, 2] = "R. albus"
pathogen_AMR <- pathogen_AMR %>%
  filter(!grepl("S. mannosilyticum|S. sp. FSL K6-2383", name))
# Seaweed
pathogen_AMR[60, 2] = "C. innocuum"
pathogen_AMR[61:63, 2] = "C. scindens"
pathogen_AMR[40, 2] = "M. formatexigens"
pathogen_AMR[32, 2] = "F. nucleatum"
pathogen_AMR[6, 2] = "C. Pelagibacter"
pathogen_AMR[21, 2] = "C. sacch."
pathogen_AMR[51, 2] = "R. intestinalis"
pathogen_AMR[64:65, 2] = "u. Subdoligranulum"

pathogen_AMR <- pathogen_AMR %>%
  filter(!grepl("S. mannosilyticum|S. sp. FSL K6-2383", name))


# Filter by species
# df <- abri_kraken2_filtered %>%
#   filter(grepl("Clostridioides difficile|Enterococcus faecium|Staphylococcus aureus|Acinetobacter baumannii|Klebsiella pneumoniae|Klebsiella quasipneumoniae|Enterobacter cloacae complex|Enterobacter hormaechei|Escherichia coli|Streptococcus pneumoniae|Clostridium acetobutylicum|Clostridium perfringens|Neisseria cinerea|Streptococcus parasanguinis|Streptococcus pseudopneumoniae|Streptococcus sp. 116-D4|Streptococcus|Haemophilus|Moraxella|Staphylococcus|Corynebacterium|Neisseria|Prevotella|Veillonella", name))
df <- pathogen_AMR

# Create a connection matrix
categories <- unique(c(sort(df$sample), sort(df$name), sort(df$GENE), sort(df$RESISTANCE)))
mat <- matrix(0, nrow = length(categories), ncol = length(categories), dimnames = list(categories, categories))

for (i in 1:nrow(df)) {
  mat[df$sample[i], df$name[i]] <- 1
  mat[df$name[i], df$GENE[i]] <- 1
  mat[df$GENE[i], df$RESISTANCE[i]] <- 1
}

# Set colors for categories
color_sectors <- c(
  paletteer_d("ggsci::default_igv", n = length(unique(df$sample))),
  paletteer_d("ggsci::default_igv", n = length(unique(df$name))),
  paletteer_d("ggsci::default_igv", n = length(unique(df$GENE))),
  paletteer_d("ggsci::default_igv", n = length(unique(df$RESISTANCE)))
)



# Graph parameters
circos.par(track.height = 0.1, 
           start.degree = 105, 
           gap.degree = 2, 
           canvas.xlim = c(-1, 1), 
           canvas.ylim = c(-1, 1), 
           circle.margin = c(1, 1), 
           unit.circle.segments = 500)

# Create the chord diagram
chordDiagram(mat, 
             transparency = 0.5, 
             annotationTrack = "grid", 
             scale = FALSE, directional = 1, 
             diffHeight = mm_h(3), 
             grid.col = color_sectors, 
             preAllocateTracks = list(track.height = 0.1, 
                                      unit.circle.segments = 105, 
                                      start.degree = 95, 
                                      scale = TRUE)
             )

# Add labels to the sectors
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.8), cex = 1)
}, bg.border = NA)



## reset the graphic parameters and internal variables
circos.clear()

# Clear environment
rm(list = ls())

# Clear packages
detach("package:datasets", unload = TRUE)

# Clear plots
dev.off()  # But only if there IS a plot

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)
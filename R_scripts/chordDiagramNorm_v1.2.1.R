# Purpose: Chord diagram linking samples -> taxa -> AMR genes -> drug-class
#          categories using the abricate+Bracken merged dataset (CARD subset).
# Input:   data/abri_kraken2Bracken_merged.csv (cols: sample, name, GENE,
#                                               RESISTANCE, DATABASE, Treatment,
#                                               sampleCount, ...)
# Output:  Chord diagram rendered to the active graphics device.

library(readr)
library(circlize)
library(dplyr)
library(tidyr)
library(stringr)
library(grid)
library(paletteer)

# Load dataset
abri_kraken2Bracken_merged <- read_csv("data/abri_kraken2Bracken_merged.csv")

# Helper: title-case and replace underscores
clean_string <- function(x) {
  str_to_title(str_replace_all(x, "_", " "))
}

# Helper: collapse multi-class RESISTANCE strings to a single label
classify_resistance <- function(res_string) {
  if (is.null(res_string) || is.na(res_string) || res_string == "") return(NA)
  classes        <- trimws(unlist(strsplit(res_string, ";")))
  classes        <- classes[classes != ""]
  unique_classes <- unique(classes)
  if (length(unique_classes) == 1) return(unique_classes)
  if (length(unique_classes) > 1) return("Multi-drug")
  return(NA)
}

# Filter and annotate CARD entries
TaxageneTable_filtered <- abri_kraken2Bracken_merged %>%
  mutate(RESISTANCE = str_replace(RESISTANCE,
                                  "lincosamide;macrolide;streptogramin;streptogramin_A;streptogramin_B",
                                  "MLS")) %>%
  mutate(RESISTANCE = sapply(RESISTANCE, clean_string)) %>%
  mutate(RESISTANCE = sapply(RESISTANCE, classify_resistance)) %>%
  mutate(RESISTANCE = str_replace(RESISTANCE, "Mls", "MLS")) %>%
  arrange(Treatment) %>%
  filter(grepl("card", DATABASE)) %>%
  select(sample, name, GENE, RESISTANCE, sampleCount, Treatment) %>%
  distinct(name, RESISTANCE, .keep_all = TRUE) %>%
  filter(!grepl("Alistipes|Azorhizobium|cellular organisms|Bacillota|Clostridia|Lachnospiraceae|Bacteria|root|Enterobacterales|Enterobacteriaceae|Terrabacteria group|Homo sapiens|Pseudomonadota|Bacteroidales|Bifidobacterium|Bacteroidota", name)) %>%
  filter(!(str_count(name, "\\S+") == 1 & name != "Enterococcus")) %>%
  filter(sampleCount > 50)

# Abbreviate species and clean gene names
df <- TaxageneTable_filtered %>%
  arrange(Treatment) %>%
  mutate(
    name = str_replace(name, "^(\\w)\\w+\\s(\\w+)", "\\1. \\2"),
    name = if_else(str_length(name) > 15,
                   str_replace(name, "^((?:\\S+\\s+){2}).*", "\\1"),
                   name)
  ) %>%
  mutate(
    name = case_when(
      grepl("\\[Clostridium\\] asparagiforme", name) ~ "C. asparagiforme",
      grepl("\\[Clostridium\\] scindens",      name) ~ "C. scindens",
      grepl("\\[Ruminococcus\\] lactaris",     name) ~ "R. lactaris ",
      grepl("\\[Ruminococcus\\] torques",      name) ~ "R. torques",
      grepl("\\[Clostridium\\] innocuum",      name) ~ "C. innocuum",
      TRUE ~ name
    ),
    name = ifelse(grepl("B. pullorum subsp. gallinarum", name), "B. pullorum", name)
  ) %>%
  mutate(
    GENE = case_when(
      grepl("vanR_gene_in_vanA_cluster", GENE) ~ "vanRA",
      grepl("vanH_gene_in_vanA_cluster", GENE) ~ "vanHA",
      grepl("vanS_gene_in_vanA_cluster", GENE) ~ "vanSA",
      TRUE ~ GENE
    ),
    name = ifelse(grepl("Rahnella aquatilis CIP 78.65 = ATCC 33071", name),
                  "Rahnella aquatilis", name)
  ) %>%
  rename(SAMPLE = sample, NAME = name) %>%
  select(SAMPLE, NAME, GENE, RESISTANCE)

# Build connection matrix
categories <- unique(c(df$SAMPLE, sort(df$NAME), sort(df$GENE), sort(df$RESISTANCE)))
mat <- matrix(0,
              nrow = length(categories),
              ncol = length(categories),
              dimnames = list(categories, categories))

for (i in seq_len(nrow(df))) {
  mat[df$SAMPLE[i], df$NAME[i]]  <- 1
  mat[df$NAME[i],   df$GENE[i]]  <- 1
  mat[df$GENE[i],   df$RESISTANCE[i]] <- 1
}

# Colours for each node layer
color_sectors <- c(
  paletteer_d("ggsci::default_ucscgb", n = length(unique(df$SAMPLE))),
  paletteer_d("ggsci::default_ucscgb", n = length(unique(df$NAME))),
  paletteer_d("ggsci::default_ucscgb", n = length(unique(df$GENE))),
  paletteer_d("ggsci::default_ucscgb", n = length(unique(df$RESISTANCE)))
)

circos.clear()
print("Plotting the chord diagram...")

circos.par(track.height = 0.1, start.degree = 152, gap.degree = 2,
           canvas.xlim = c(-1, 1), canvas.ylim = c(-1, 1),
           circle.margin = c(1, 1), unit.circle.segments = 500)

chordDiagram(mat, transparency = 0.5, annotationTrack = "grid",
             scale = FALSE, directional = 1, diffHeight = mm_h(3),
             grid.col = color_sectors,
             preAllocateTracks = list(track.height = 0.1,
                                      unit.circle.segments = 100,
                                      start.degree = 90,
                                      scale = TRUE))

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector.index,
              facing = "clockwise", niceFacing = TRUE,
              adj = c(0, 0.8), cex = 1)
}, bg.border = NA)

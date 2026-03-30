# Purpose: Venn diagram showing species shared across treatment groups (Reference
#          diet, Soyabean meal, Dulce) based on Bracken taxonomic counts.
# Input:   data/bracken_arranged.csv  (cols: #perc, tot_all, tot_lvl, name, taxid,
#                                      D19..D36 [sample columns])
#          data/chicken_metadata.csv  (cols: sample, Treatment, ...)
# Output:  Venn diagram rendered to the active graphics device.

library(VennDiagram)
library(dplyr)
library(readr)
library(tidyr)
library(tidyverse)
library(stringr)

# Load datasets
bracken_arranged  <- read_csv("data/bracken_arranged.csv")
chicken_metadata  <- read_csv("data/chicken_metadata.csv")
chicken_metadata  <- na.omit(chicken_metadata)

# Standardise treatment labels
modify_strings <- function(df, column_name, target_string, replacement_string) {
  df %>% mutate(!!sym(column_name) := str_replace_all(
    !!sym(column_name), fixed(target_string), replacement_string))
}

chicken_metadata <- modify_strings(chicken_metadata, "Treatment",
                                   "Reference diet (Wheat -soyabean)", "Reference diet")
chicken_metadata <- modify_strings(chicken_metadata, "Treatment",
                                   "Seaweed", "Dulce")

# Pivot to long format, filter contaminants and rare taxa
bracken_pivoted <- bracken_arranged %>%
  pivot_longer(cols = 4:21, names_to = "sample", values_to = "count") %>%
  filter(!grepl("root|Homo sapiens|cellular organisms|unclassified|Bacteria|environmental samples", name)) %>%
  filter(!is.na(name)) %>%
  group_by(sample) %>%
  filter(count > 5) %>%
  mutate(Percentage = (count / sum(count)) * 100) %>%
  ungroup()

bracken_merged <- merge(bracken_pivoted, chicken_metadata, by = "sample")

# Build treatment x species count matrix (species-level only, abbreviated names)
taxaTreatment_dataset <- bracken_merged %>%
  filter(!grepl(
    "Terrabacteria group|Bacteroidota/Chlorobiota group|FCB group|Cyanobacteriota/Melainabacteria group|Rhizobium/Agrobacterium group|Candida/Lodderomyces clade|Sinorhizobium/Ensifer group|Plasmodium",
    name
  )) %>%
  filter(str_count(name, "\\S+") > 1) %>%
  mutate(
    name = str_replace(name, "^(\\w)\\w+\\s(\\w+)", "\\1. \\2"),
    name = if_else(
      str_length(name) > 30,
      str_replace(name, "^((?:\\S+\\s+){2}).*", "\\1"),
      name
    )
  ) %>%
  mutate(
    name = case_when(
      grepl("\\[Clostridium\\] asparagiforme ", name) ~ "C. asparagiforme",
      grepl("\\[Clostridium\\] scindens", name)       ~ "C. scindens",
      grepl("\\[Ruminococcus\\] lactaris", name)      ~ "R. lactaris ",
      grepl("\\[Ruminococcus\\] torques", name)       ~ "R. torques",
      grepl("\\[Clostridium\\] innocuum", name)       ~ "C. innocuum",
      grepl("\\[Clostridium\\] colinum", name)        ~ "C. colinum",
      grepl("\\[Clostridium\\] hylemonae", name)      ~ "C. hylemonae",
      TRUE ~ name
    ),
    name = ifelse(grepl("u. Subdoligranulum", name), "U. Subdoligranulum", name)
  ) %>%
  select(name, Treatment, count) %>%
  tidyr::pivot_wider(
    names_from  = Treatment,
    values_from = count,
    values_fill = 0,
    values_fn   = sum
  ) %>%
  distinct(name, .keep_all = TRUE) %>%
  column_to_rownames("name")

# Convert to binary matrix and create Venn sets
treatment_dataset_matrix        <- data.matrix(taxaTreatment_dataset, rownames.force = NA)
treatment_dataset_matrix_Bmatrix <- as.matrix((treatment_dataset_matrix > 0) + 0)

sets       <- apply(treatment_dataset_matrix_Bmatrix, 2, function(col) which(col == 1))
names(sets) <- colnames(treatment_dataset_matrix_Bmatrix)
print(sets)

# Venn diagram
venn.plot <- venn.diagram(
  x              = sets,
  category.names = c("Reference diet", "Soyabean meal", "Dulce"),
  fill           = c(`Reference diet` = "#c2320e", `Soyabean meal` = "#c2980e", `Dulce` = "#04540a"),
  alpha          = 0.7,
  height         = 50,
  width          = 50,
  filename       = NULL,
  cat.cex        = 1.3,
  cex            = 2
)
grid.newpage()
grid.draw(venn.plot)

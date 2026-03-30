# Purpose: Venn diagram of taxon IDs (by DATABASE category) shared across the
#          three treatment groups using the abricate+Bracken merged dataset.
# Input:   data/abri_kraken2Bracken_merged.csv (cols: sample, name, taxid, GENE,
#                                               DATABASE, Treatment, ...)
# Output:  Venn diagram rendered to the active graphics device.

library(VennDiagram)
library(dplyr)
library(readr)
library(tidyr)
library(tidyverse)
library(stringr)

# Load dataset
abri_kraken2_merged <- read_csv("data/abri_kraken2Bracken_merged.csv")

# Filter to AMR/VF/MGE annotations, species-level taxa only
abri_kraken2_merged <- abri_kraken2_merged %>%
  filter(grepl("card|vfdb|plasmidfinder", DATABASE)) %>%
  filter(!grepl("Terrabacteria group|Bacteroidota/Chlorobiota group|FCB group", name)) %>%
  filter(str_count(name, "\\S+") > 1)

# Split by treatment group, keeping unique taxon IDs
Wheat_soyabean_dataset <- abri_kraken2_merged %>%
  filter(Treatment == "Reference diet") %>%
  distinct() %>%
  rename(Wheat_soyabean = Treatment) %>%
  select(taxid, Wheat_soyabean)

Seaweed_dataset <- abri_kraken2_merged %>%
  filter(Treatment == "Dulce") %>%
  distinct() %>%
  rename(Seaweed = Treatment) %>%
  select(taxid, Seaweed)

Soyabean_meal_dataset <- abri_kraken2_merged %>%
  filter(Treatment == "Soyabean meal") %>%
  distinct() %>%
  rename(Soyabean_meal = Treatment) %>%
  select(taxid, Soyabean_meal)

Wheat_soyabean_unique  <- length(unique(Wheat_soyabean_dataset$taxid))
Seaweed_unique         <- length(unique(Seaweed_dataset$taxid))
Soyabean_meal_unique   <- length(unique(Soyabean_meal_dataset$taxid))

# Merge all treatment datasets
treatment_dataset  <- merge(Wheat_soyabean_dataset, Seaweed_dataset, by = "taxid", all = TRUE)
treatment_dataset1 <- merge(treatment_dataset, Soyabean_meal_dataset, by = "taxid", all = TRUE)

treatment_dataset1 <- treatment_dataset1 %>%
  distinct(taxid, .keep_all = TRUE) %>%
  remove_rownames() %>%
  column_to_rownames(var = "taxid") %>%
  mutate(across(everything(), ~ if_else(. == "", NA_character_, .))) %>%
  filter(!if_all(everything(), is.na))

# Convert to binary matrix and create Venn sets
treatment_dataset_matrix        <- data.matrix(treatment_dataset1, rownames.force = NA) %>%
  replace(is.na(.), 0)
sum(is.na(treatment_dataset_matrix))

treatment_dataset_matrix_Bmatrix <- as.matrix((treatment_dataset_matrix > 0) + 0)

sets       <- apply(treatment_dataset_matrix_Bmatrix, 2, function(col) which(col == 1))
names(sets) <- colnames(treatment_dataset_matrix_Bmatrix)
print(sets)

# Venn diagram
venn.plot <- venn.diagram(
  x              = sets,
  category.names = c("Reference diet", "Dulce", "Soyabean meal"),
  fill           = c(Wheat_soyabean = "#c2320e", Seaweed = "#04540a", Soyabean_meal = "#c2980e"),
  alpha          = 0.3,
  height         = 50,
  width          = 50,
  filename       = NULL,
  cat.cex        = 1.2,
  cex            = 2
)
grid.newpage()
grid.draw(venn.plot)

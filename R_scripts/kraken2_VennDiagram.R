# Purpose: Venn diagram of taxon IDs shared across the three treatment groups
#          (Reference diet, Seaweed, Soyabean meal) using the Kraken2 filtered
#          dataset. Identifies core, shared, and treatment-unique taxa.
# Input:   ../data/abri_kraken2_filtered.csv (cols: SAMPLE, TAXID, NAME, GENE,
#                                             DATABASE, TREATMENT, ...)
# Output:  Venn diagram rendered to the active graphics device.

library(VennDiagram)
library(dplyr)
library(readr)
library(tidyr)
library(tidyverse)
library(stringr)

# Load dataset
abri_kraken2_filtered <- read_csv("../data/abri_kraken2_filtered.csv")

# Split by treatment group, keeping unique taxon IDs
Reference_dataset <- abri_kraken2_filtered %>%
  filter(TREATMENT == "Reference diet") %>%
  distinct() %>%
  rename(Reference_diet = TREATMENT) %>%
  select(TAXID, Reference_diet)

Seaweed_dataset <- abri_kraken2_filtered %>%
  filter(TREATMENT == "Seaweed") %>%
  distinct() %>%
  rename(Seaweed = TREATMENT) %>%
  select(TAXID, Seaweed)

Soyabean_dataset <- abri_kraken2_filtered %>%
  filter(TREATMENT == "Soyabean meal") %>%
  distinct() %>%
  rename(Soyabean_meal = TREATMENT) %>%
  select(TAXID, Soyabean_meal)

Reference_unique  <- length(unique(Reference_dataset$TAXID))
Seaweed_unique    <- length(unique(Seaweed_dataset$TAXID))
Soyabean_unique   <- length(unique(Soyabean_dataset$TAXID))

# Merge all treatment datasets
treatment_dataset  <- merge(Reference_dataset, Seaweed_dataset, by = "TAXID", all = TRUE)
treatment_dataset1 <- merge(treatment_dataset, Soyabean_dataset, by = "TAXID", all = TRUE)

treatment_dataset1 <- treatment_dataset1 %>%
  distinct(TAXID, .keep_all = TRUE) %>%
  remove_rownames() %>%
  column_to_rownames(var = "TAXID") %>%
  mutate(across(everything(), ~ if_else(. == "", NA_character_, .))) %>%
  filter(!if_all(everything(), is.na))

# Convert to binary matrix and create Venn sets
treatment_dataset_matrix <- data.matrix(treatment_dataset1, rownames.force = NA) %>%
  replace(is.na(.), 0)
sum(is.na(treatment_dataset_matrix))

treatment_dataset_matrix_Bmatrix <- as.matrix((treatment_dataset_matrix > 0) + 0)

sets        <- apply(treatment_dataset_matrix_Bmatrix, 2, function(col) which(col == 1))
names(sets) <- colnames(treatment_dataset_matrix_Bmatrix)
print(sets)

# Venn diagram
venn.plot <- venn.diagram(
  x              = sets,
  category.names = c("Reference diet", "Seaweed", "Soyabean meal"),
  fill           = c(Reference_diet = "#c2320e", Seaweed = "#04540a", Soyabean_meal = "#c2980e"),
  alpha          = 0.3,
  height         = 50,
  width          = 50,
  filename       = NULL,
  cat.cex        = 1.2,
  cex            = 2
)
grid.newpage()
grid.draw(venn.plot)

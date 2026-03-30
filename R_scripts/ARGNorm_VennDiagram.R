# Purpose: Venn diagram of ARG drug-class overlap across the three treatment groups
#          using TPM-normalised gene-table data.
# Input:   data/genetable_normdata.csv (cols: sample, GENE, PRODUCT, RESISTANCE,
#                                       DATABASE, TPM, Treatment, ...)
# Output:  Venn diagram rendered to the active graphics device.

library(VennDiagram)
library(dplyr)
library(readr)
library(tidyr)
library(tidyverse)
library(stringr)

# Load dataset
genetable_normdata <- read_csv("data/genetable_normdata.csv")

# Helper functions
clean_string <- function(x) str_to_title(str_replace_all(x, "_", " "))

classify_resistance <- function(res_string) {
  if (is.null(res_string) || is.na(res_string) || res_string == "") return(NA)
  classes        <- trimws(unlist(strsplit(res_string, ";")))
  classes        <- classes[classes != ""]
  unique_classes <- unique(classes)
  if (length(unique_classes) == 1) return(unique_classes)
  if (length(unique_classes) > 1) return("Multi-drug")
  return(NA)
}

# Filter to ARGs; classify resistance; pivot to treatment x drug-class matrix
geneTable_filtered <- genetable_normdata %>%
  mutate(RESISTANCE = str_replace(
    RESISTANCE,
    "lincosamide;macrolide;streptogramin;streptogramin_A;streptogramin_B", "MLS"
  )) %>%
  mutate(RESISTANCE = sapply(RESISTANCE, clean_string)) %>%
  mutate(RESISTANCE = sapply(RESISTANCE, classify_resistance)) %>%
  mutate(RESISTANCE = str_replace(RESISTANCE, "Mls", "MLS")) %>%
  arrange(Treatment) %>%
  filter(grepl("card", DATABASE)) %>%
  mutate(TPM_log_count = log(TPM + 1)) %>%
  select(RESISTANCE, Treatment, TPM_log_count) %>%
  tidyr::pivot_wider(
    names_from  = Treatment,
    values_from = TPM_log_count,
    values_fill = 0,
    values_fn   = sum
  ) %>%
  distinct(RESISTANCE, .keep_all = TRUE) %>%
  column_to_rownames("RESISTANCE")

# Convert to binary matrix and create Venn sets
treatment_dataset_matrix <- data.matrix(geneTable_filtered, rownames.force = NA)

sum(is.na(treatment_dataset_matrix))

treatment_dataset_matrix_Bmatrix <- as.matrix((treatment_dataset_matrix > 0) + 0)

sets       <- apply(treatment_dataset_matrix_Bmatrix, 2, function(col) which(col == 1))
names(sets) <- colnames(treatment_dataset_matrix_Bmatrix)
print(sets)

# Venn diagram
venn.plot <- venn.diagram(
  x              = sets,
  category.names = c("Reference diet", "Dulce", "Soyabean meal"),
  fill           = c(Reference_diet = "#c2320e", Dulce = "#04540a", Soyabean_meal = "#c2980e"),
  alpha          = 0.6,
  height         = 50,
  width          = 50,
  filename       = NULL,
  cat.cex        = 1.5,
  cex            = 2
)
grid.newpage()
grid.draw(venn.plot)

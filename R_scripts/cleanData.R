
# Load libraries
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(purrr)


# Clean Abricate summary report
summary_reportClean <- read_csv("../data/summary_reportClean.csv")
# View(summary_reportClean)


# Kraken combined report
kraken2_db1_combined_reports <- read_delim("../data/kraken2_db1_combined_reports.txt",
                                               delim = "\t", escape_double = FALSE,
                                               trim_ws = TRUE)

# Bracken combined report
Bracken_db2_combined_reports <- read_delim("../data/kraken2_db2-bracken_combined_reports.txt",
                                           delim = "\t", escape_double = FALSE,
                                           trim_ws = TRUE)

# Recentrifuge contaminats list
contaminats_data <- read_csv("../data/recentrifuge_contaminants.csv")



# Subset dataset
taxaCount_Kraken2 <- kraken2_db1_combined_reports %>%
  select("#perc", "taxid", "name", "tot_all")

# Rename the column to "tot_all" to "count"
names(taxaCount_Kraken2)[names(taxaCount_Kraken2) == "tot_all"] <- "count"

# taxaCount_Bracken <- Bracken_db2_combined_reports %>%
#   select("#perc", "name")


# Store Abricate summary cleaned into summary dataframe
Summary <- as.data.frame(summary_reportClean)

# Rename the column to "sample"
names(Summary)[names(Summary) == "#FILE"] <- "sample"


# Remove the ".fasta" extension
Summary$sample <- gsub(".fasta", "", Summary$sample)

# Remove everything after the first underscore
Summary$sample <- sub("_.*", "", Summary$sample)

# Step 1: Split the string into two parts at the underscore
Summary <- Summary %>%
  separate(SEQUENCE, into = c("sequence", "taxid"), sep = "_")

# Step 2: Extract the part after the '|' and create a new column called 'taxid'
Summary <- Summary %>%
  mutate(taxid = sub(".*\\|", "", taxid))

# Sample category
sample_category <- as.factor(Summary$sample)
taxa_id <- as.factor(Summary$taxid)

# Merge Abricate and Kraken2 output
abri_kraken2_merged <- merge(Summary,taxaCount_Kraken2, by = c("taxid"))

# Clean-up the abri_kraken2_merged dataset
abri_kraken2_cleaned <- abri_kraken2_merged %>%
  filter(!grepl("root|Homo sapiens|cellular organisms|unclassified|Bacteria|environmental samples", name))

# Save abri_kraken2_merged dadaset as csv file
write.csv(abri_kraken2_cleaned, "../data/abri_kraken2_filtered.csv", row.names = FALSE)


# Bracken report clean-up
colnames(Bracken_db2_combined_reports)
bracken_report <- Bracken_db2_combined_reports %>%
  select("#perc",   "tot_all",  "tot_lvl",  "D25", "D32", "D22", "D35", "D28", 
         "D21", "D29", "D31", "D27", "D20", "D34", "D24", "D36", "D23", "D33", 
         "D19", "D26", "D30", "taxid", "name")


# Save bracken report arranged
write.csv(bracken_report, "../data/bracken_arranged.csv", row.names = FALSE)



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

# Load the libraries
library(readr)
library(dplyr)
library(tidyverse)
library(reshape2)
library(stringr)

# Load dataset
abri_kraken2Bracken_merged <- read_csv("../data/abri_kraken2Bracken_merged.csv")

# TPM Calculation Steps

# Calculate Reads Per Kilobase (RPK): Divide the raw read count for each gene by its length in kilobases (kb).
df_genetable <- abri_kraken2Bracken_merged %>%
  group_by(sample, GENE) %>%
  mutate(GeneLength = abs(END - START) + 1) %>%
  mutate(geneCount = n()) %>% 
  summarise(
    meanLength = round(mean(GeneLength), 1),
    meanCount = round(mean(geneCount), 1),
    .groups = "drop"
  ) %>%
  mutate(RPK = meanCount/(meanLength/1000)) %>% # Extract and Calculate RPK (Reads Per Kilobase)
  select(sample, GENE, meanLength, meanCount, RPK) %>%
  distinct()

# Calculate Scaling Factor: Sum all the RPK values in the sample and divide by 1,000,000.
# Divide each gene's RPK by the scaling factor. 
df_genetable_scaled <- df_genetable %>%
  group_by(sample) %>%
  mutate(TPM = (RPK/(sum(RPK)/ 1e6))) %>%
  ungroup()

# Double check TPM scaling
check_TPM <- df_genetable_scaled %>%
  group_by(sample) %>%
  #filter(sample == "D19") %>%
  mutate(totalTPM = sum(TPM)) %>%
  mutate(totalRPK = sum(RPK)) %>%
  mutate(ScaFactor = totalRPK/totalTPM) %>%
  select(sample, totalRPK, totalTPM, ScaFactor) %>%
  distinct()


# Add metadata to Gene table scaled
Metadata = abri_kraken2Bracken_merged %>%
  select(sample, GENE, PRODUCT, RESISTANCE, DATABASE, Treatment)%>%
  distinct(sample, GENE, .keep_all = T)

# Merge Gene table scaled and metadata
df_genetable_scaledMeta <- merge(df_genetable_scaled, Metadata, by=c("sample", "GENE"))

# Number of Unique Genes
uniqueGene <- length(unique(df_genetable_scaledMeta$GENE))

# Save Gene norm table 
write.csv(df_genetable_scaledMeta, "../data/genetable_normdata.csv", row.names = F)


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


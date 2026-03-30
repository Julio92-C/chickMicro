# Purpose: Summary statistics and bar/stacked plots for AMR, VFs, and MGEs
#          detected across treatment groups in broiler chicken metagenomics.
# Input:   data/abri_kraken2_cleaned.csv  (cols: sample, name, GENE, RESISTANCE,
#                                          DATABASE, sequence, START, END, ...)
#          data/chicken_metadata.csv      (cols: sample, Treatment, ...)
# Output:  Plots rendered to the active graphics device (no files written).

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(paletteer)
library(stringr)

# Load datasets
abri_kraken2_filtered <- read_csv("data/abri_kraken2_cleaned.csv")
chicken_metadata      <- read_csv("data/chicken_metadata.csv")
chicken_metadata      <- na.omit(chicken_metadata)

# Standardise treatment label
modify_strings <- function(df, column_name, target_string, replacement_string) {
  df %>% mutate(!!sym(column_name) := str_replace_all(
    !!sym(column_name), fixed(target_string), replacement_string))
}

chicken_metadata <- modify_strings(
  chicken_metadata, "Treatment",
  "Reference diet (Wheat -soyabean)", "Reference diet"
)

# Merge abricate output with metadata
abri_kraken2_merged <- merge(abri_kraken2_filtered, chicken_metadata, by = "sample")

# Remove high-level / human taxonomy rows
abri_kraken2_filtered <- abri_kraken2_merged %>%
  filter(name != "root") %>%
  filter(name != "Bacteria") %>%
  filter(name != "Homo sapiens")

# Summary counts
sample     <- length(unique(abri_kraken2_filtered$sample))
taxa       <- length(unique(abri_kraken2_filtered$name))
gene       <- length(unique(abri_kraken2_filtered$GENE))
resistance <- length(unique(abri_kraken2_filtered$RESISTANCE))

# Taxa with >5 occurrences
df_species <- abri_kraken2_merged %>%
  group_by(name) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  filter(count > 5) %>%
  arrange(desc(count))

Species_count <- ggplot(df_species, aes(x = reorder(name, -count), y = count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = round(count / 1, 1)), vjust = -1.5, size = 2.8) +
  labs(title = "", x = "Taxonomy Level", y = "log10(Taxa Count)") +
  scale_y_log10() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14))

Species_count

# Relative abundance by treatment
species_treatment <- abri_kraken2_filtered %>%
  group_by(Treatment, name) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(Percentage = (count / sum(count)) * 100)

taxa_f   <- as.factor(species_treatment$name)
t        <- length(levels(taxa_f))
t_palette <- paletteer_d("palettesForR::Visibone", n = t)

ggplot(species_treatment, aes(x = factor(Treatment), y = Percentage, fill = name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = t_palette) +
  labs(x = "Treatment Groups", y = "Percentage", fill = "Taxon") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 16))

# Deduplicated species-level subset
abri_kraken2_uniqued <- abri_kraken2_filtered %>%
  filter(!grepl("Terrabacteria group|Bacteroidota/Chlorobiota group|FCB group", name)) %>%
  filter(str_count(name, "\\S+") > 1) %>%
  distinct(GENE, sequence, .keep_all = TRUE)

uniSeq    <- length(unique(abri_kraken2_uniqued$sequence))
uniTaxa   <- length(unique(abri_kraken2_uniqued$name))
uniGene   <- length(unique(abri_kraken2_uniqued$GENE))
uniSample <- length(unique(abri_kraken2_uniqued$sample))

# AMR (CARD) gene counts
df_AMR <- abri_kraken2_uniqued %>%
  filter(DATABASE == "card") %>%
  group_by(GENE) %>%
  summarise(Genecount = n(), .groups = "drop") %>%
  filter(Genecount > 2) %>%
  arrange(desc(Genecount))

df_AMR[10, 1] <- "rpoB_mutants"
df_AMR[12, 1] <- "ileS"

AMR_count <- ggplot(df_AMR, aes(x = reorder(GENE, -Genecount), y = Genecount)) +
  geom_bar(stat = "identity", fill = "#b02d02") +
  geom_text(aes(label = round(Genecount / 1, 1)), vjust = -0.5, size = 2.8) +
  labs(title = "", x = "Antibiotics Resistance Genes (ARGs)", y = "ARG Count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 16))

AMR_count

# VFs (VFDB) — extract functional category from product description
extract_productFunction <- function(string) {
  match  <- regmatches(string, regexpr("\\- ([^\\(]+) \\(", string))
  result <- gsub("\\- | \\(", "", match)
  return(result)
}

refactored_df <- abri_kraken2_filtered %>%
  filter(DATABASE == "vfdb") %>%
  mutate(Functions = sapply(PRODUCT, extract_productFunction))

unique_seqVFs <- length(unique(refactored_df$sequence))

vfdb_filtered <- refactored_df %>%
  group_by(Functions, Treatment) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(count > 1) %>%
  arrange(desc(count))

uniqueFunction <- length(unique(vfdb_filtered$Functions))

Functions_f <- as.factor(vfdb_filtered$Functions)
f           <- length(levels(Functions_f))
f_palette   <- paletteer_d("ggsci::default_nejm", n = f)

VFs_count <- ggplot(vfdb_filtered, aes(x = reorder(Functions, -count), y = count)) +
  geom_bar(stat = "identity", aes(fill = Treatment)) +
  scale_fill_manual(values = t_palette) +
  labs(title = "", x = "Virulence Factors (VFs)", y = "VFs Count") +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 16)) +
  facet_wrap(~ Functions, scales = "free_x", nrow = 3)

VFs_count

# MGEs (PlasmidFinder) gene counts by treatment
plasmid_filtered <- abri_kraken2_filtered %>%
  filter(DATABASE == "plasmidfinder") %>%
  group_by(GENE, Treatment) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(count > 5) %>%
  arrange(desc(count))

taxa_p    <- as.factor(plasmid_filtered$GENE)
t         <- length(levels(taxa_p))
t_palette <- paletteer_d("ggsci::default_nejm", n = t)

MGEs_count <- ggplot(plasmid_filtered, aes(x = reorder(Treatment, -count), y = count)) +
  geom_bar(stat = "identity", aes(fill = GENE)) +
  scale_fill_manual(values = t_palette) +
  labs(title = "", x = "Treatment groups", y = "log10(MGE Count)") +
  scale_y_log10() +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        text = element_text(size = 15))

MGEs_count

# Database composition per sample
data_summary <- abri_kraken2_filtered %>%
  filter(grepl("card|vfdb|plasmidfinder", DATABASE)) %>%
  group_by(sample, DATABASE) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(Percentage = (count / sum(count)) * 100)

ggplot(data_summary, aes(x = factor(sample), y = Percentage, fill = DATABASE)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#b02d02", "#5c0404", "#c47a02")) +
  labs(x = "Sample Groups", y = "Percentage", fill = "DATABASE") +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 16))

# Gene count by database (pie chart)
Gene_count <- aggregate(count ~ DATABASE, data = data_summary, sum)

pie(Gene_count$count,
    labels = paste0(Gene_count$DATABASE, " / ", Gene_count$count, " / ",
                    round(100 * Gene_count$count / sum(Gene_count$count), 2), "%"),
    main   = "Genetic elements count",
    col    = c("#b02d02", "#5c0404", "#c47a02"),
    cex    = 1,
    radius = 1,
    xlim   = c(-1.5, 1.5))

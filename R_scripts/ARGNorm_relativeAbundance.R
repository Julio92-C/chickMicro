# Purpose: Stacked relative-abundance bar plots and ranked bar charts of ARG
#          drug classes derived from TPM-normalised data.
# Input:   data/genetable_normdata.csv (cols: sample, GENE, PRODUCT, RESISTANCE,
#                                       DATABASE, TPM, Treatment, ...)
# Output:  Plots rendered to the active graphics device.

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(purrr)
library(paletteer)
library(plotly)
library(htmlwidgets)

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

# Filter to ARGs; standardise MLS label; classify drug resistance
genetable_filtered <- genetable_normdata %>%
  filter(DATABASE == "card") %>%
  select(sample, GENE, PRODUCT, RESISTANCE, TPM, Treatment)

genetable_arranged <- genetable_filtered %>%
  mutate(RESISTANCE = str_replace(
    RESISTANCE,
    "lincosamide;macrolide;streptogramin;streptogramin_A;streptogramin_B", "MLS"
  ))

Resistance_category <- genetable_arranged %>%
  distinct(sample, GENE, RESISTANCE, .keep_all = TRUE) %>%
  mutate(RESISTANCE = sapply(RESISTANCE, clean_string)) %>%
  mutate(RESISTANCE = sapply(RESISTANCE, classify_resistance)) %>%
  group_by(sample) %>%
  mutate(Percentage = (TPM / sum(TPM)) * 100) %>%
  ungroup() %>%
  mutate(RESISTANCE = str_replace(RESISTANCE, "Mls", "MLS"))

Total_Perc <- sum(Resistance_category$Percentage)

resistCategory <- as.factor(Resistance_category$RESISTANCE)
r              <- length(levels(resistCategory))
r_palette      <- paletteer_d("ggsci::default_ucscgb", n = r)

# Stacked relative-abundance bar plot
ra <- ggplot(Resistance_category, aes(x = factor(sample), y = Percentage, fill = RESISTANCE)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = r_palette) +
  labs(x = "Sample Groups by Treatment", y = "Percentage (Relative abundance)") +
  theme_classic() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14)) +
  guides(fill = guide_legend(nrow = 10)) +
  labs(fill = "Drug classes") +
  facet_wrap(~ Treatment, scales = "free_x", nrow = 1)

ra

# Ranked total TPM per drug class
card_filtered <- Resistance_category %>%
  group_by(RESISTANCE) %>%
  summarise(Total_TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total_TPM))

Drug_count <- ggplot(card_filtered, aes(x = reorder(RESISTANCE, Total_TPM),
                                         y = Total_TPM, fill = RESISTANCE)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = median(card_filtered$Total_TPM), linetype = "dashed", color = "black") +
  geom_text(aes(label = paste0(round(Total_TPM / 1000, 1), "K")),
            vjust = 0.5, hjust = 1.0, size = 3.5) +
  scale_fill_manual(values = r_palette) +
  scale_y_log10() +
  labs(title = "", x = "Drug Classes", y = "Total Count") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        text = element_text(size = 14)) +
  coord_flip()

Drug_count

# Ranked total TPM per gene (drug-class on x-axis)
Resistance_category_aggregated <- genetable_arranged %>%
  distinct(sample, GENE, RESISTANCE, .keep_all = TRUE) %>%
  mutate(RESISTANCE = sapply(RESISTANCE, clean_string)) %>%
  mutate(RESISTANCE = sapply(RESISTANCE, classify_resistance)) %>%
  group_by(GENE) %>%
  summarise(Total_TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total_TPM))

AMR_count <- ggplot(Resistance_category_aggregated,
                    aes(x = reorder(GENE, -Total_TPM), y = Total_TPM)) +
  geom_bar(stat = "identity", fill = "#b02d02") +
  geom_text(aes(label = round(Total_TPM / 10000, 1)), vjust = -0.5, size = 2.8) +
  labs(title = "", x = "Drug classes", y = "log10(Total Count (TPM))") +
  scale_y_log10() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 16))

AMR_count

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyverse)
library(stringr)
library(ggsignif)  # For adding statistical significance
library(multcompView)



## Load clean dataset
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

# Merge abri_kraken2_clean and metadata dataframe
abri_kraken2_merged <- merge(abri_kraken2_clean, chicken_metadata, by ="sample")


# Filter the data set by AMR, VFs and MGEs
bracken_filtered <- abri_kraken2_merged %>%
  arrange(Treatment) %>%
  group_by(Treatment, name) %>%
  summarise(Taxa_count = n(), .groups = 'drop') %>%
  ungroup() %>% 
  mutate(Taxa_log_count = log(Taxa_count + 1))  # Adding 1 to avoid log(0)


# one-way ANOVA analysis
anova_result <- aov(Taxa_log_count ~ Treatment, data = bracken_filtered)
summary(anova_result)

# Tukey HSD post-hoc test
# tukey_result <- TukeyHSD(anova_result)
# print(tukey_result)
# 
# 
# # Extract compact letter display for group significance
# tukey_letters <- multcompLetters4(anova_result, tukey_result)
# letters_df <- as.data.frame.list(tukey_letters$Treatment)
# letters_df$Treatment <- rownames(letters_df)



# Merge letters with original data
# bracken_filtered <- left_join(bracken_filtered, letters_df, by = "Treatment")


# Generate the violin plot with statistical significance for taxa
ggplot(bracken_filtered, aes(x = Treatment, y = Taxa_log_count, fill = Treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                               "Seaweed" = "#2ca02c",
                               "Soyabean meal" = "#ff7f0e")) +
  labs(x = "Sample groups", y = "Log(Number of Taxa)") +
  theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16)) +
  #geom_text(aes(label = Letters), vjust = -0.5, size = 5) +
  annotate("text", x = 2, y = max(bracken_filtered$Taxa_log_count) + 1.3,
           label = paste("ANOVA p =", format(summary(anova_result)[[1]][["Pr(>F)"]][1], digits = 2)),
           size = 5, color = "black")


#### Filter the data set by VFs #####
abri_kraken2_filtered <- abri_kraken2_merged %>%
  filter(grepl("vfdb", DATABASE)) %>%
  arrange(name) %>%
  group_by(Treatment, GENE) %>%
  summarise(Gene_count = n(), .groups = 'drop') %>%
  ungroup() %>% 
  mutate(Gene_log_count = log(Gene_count + 1))  # Adding 1 to avoid log(0)

# one-way ANOVA analysis
anova_result <- aov(Gene_log_count ~ Treatment, data = abri_kraken2_filtered)
summary(anova_result)


# Generate the violin plot with statistical significance for VFs
ggplot(abri_kraken2_filtered, aes(x = Treatment, y = Gene_count, fill = Treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                               "Seaweed" = "#2ca02c",
                               "Soyabean meal" = "#ff7f0e")) +
  labs(x = "Sample groups", y = "Number of VFs") +
  theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16)) +
  #geom_text(aes(label = Letters), vjust = -0.5, size = 5) +
  annotate("text", x = 3, y = max(abri_kraken2_filtered$Gene_count) + 1.5,
           label = paste("ANOVA p =", format(summary(anova_result)[[1]][["Pr(>F)"]][1], digits = 2)),
           size = 5, color = "black")




#### Filter the data set by MGEs #####
abri_kraken2_filtered <- abri_kraken2_merged %>%
  filter(grepl("plasmidfinder", DATABASE)) %>%
  arrange(name) %>%
  group_by(Treatment, GENE) %>%
  summarise(Gene_count = n(), .groups = 'drop') %>%
  ungroup() %>% 
  mutate(Gene_log_count = log(Gene_count + 1))  # Adding 1 to avoid log(0)

# one-way ANOVA analysis
anova_result <- aov(Gene_log_count ~ Treatment, data = abri_kraken2_filtered)
summary(anova_result)


# Generate the violin plot with statistical significance for AMR
ggplot(abri_kraken2_filtered, aes(x = Treatment, y = Gene_log_count, fill = Treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("Reference diet" = "#BC3C29FF",
                               "Seaweed" = "#0072B5FF",
                               "Soyabean meal" = "#E18727FF")) +
  labs(x = "Sample groups", y = "Log(Number of MGEs)") +
  theme_classic() +
  theme(legend.position = "none", text = element_text(size = 16)) +
  #geom_text(aes(label = Letters), vjust = -0.5, size = 5) +
  annotate("text", x = 3, y = max(abri_kraken2_filtered$Gene_log_count) + 1.8,
           label = paste("ANOVA p =", format(summary(anova_result)[[1]][["Pr(>F)"]][1], digits = 2)),
           size = 5, color = "black")



# CLEAN UP #################################################

# Clear environment
rm(list = ls()) 

# Clear plots
if (length(dev.list()) > 0) {
  dev.off()  # Only if there is an active plot
}

# Clear console
cat("\014")  # ctrl+L






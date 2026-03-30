# Load libraries
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(paletteer)
library(gtsummary)
library(gt)
library(glue)



# Load dataset
abri_kraken2_filtered <- read_csv("../data/abri_kraken2_filtered.csv")

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


# Merge dataset and metadata dataframe
abri_kraken2_merged <- merge(abri_kraken2_filtered, chicken_metadata, by = "sample")

# Perform one-way ANOVA to assess statistical differences between treatment groups
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


# Extract and format the p-value
anova_p <- format(summary(anova_result)[[1]][["Pr(>F)"]][1], digits = 2)

# Summarise gene counts per treatment
data_stats = (abri_kraken2_merged %>% 
                filter(DATABASE == "plasmidfinder") %>%
                rename(
                  "Coverage" = `%COVERAGE`,
                  "Identity"  = `%IDENTITY`
                       ) %>% 
                # build a gtsummary table
                tbl_summary(
                  by = Treatment,
                  include = c("GENE"),
                  statistic = list(
                    all_continuous() ~ "{mean} ({sd})",
                    all_categorical() ~ "{n} / {N} ({p}%)"
                  ),
                  digits = all_continuous() ~ 2,
                  label = c(GENE ~ "MGEs"),
                  sort = list(everything() ~ "frequency"),
                  missing_text = "Missing",
                  missing = "no"
                )) %>%
  #add_p() %>%
  add_overall() %>%
  # CONVERT TO A {gt} TABLE! VERY IMPORTANT STEP!
  as_gt() %>%
  tab_header(md(glue("**Table 1. Metagenomics Statistics Results By Treatment (_ANOVA p = {anova_p}_)**")))


# Print the summary table
print(data_stats)

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

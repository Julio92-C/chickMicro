
# Load libraries
library(vegan)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(plotly)
library(paletteer)
library(funrar)
library(gtsummary)
library(gt)
library(glue)


# Load dataset
genetable_normdata <- read_csv("../data/genetable_normdata.csv")

# Metadata 
metadata <- genetable_normdata %>%
  distinct(sample, .keep_all = T) %>%
  select(sample, Treatment)


#### Filter the data set by MGEs #####
diversity_filtered <- genetable_normdata %>%
  #filter(grepl("plasmidfinder", DATABASE)) %>%
  #filter(grepl("card", DATABASE)) %>%
  filter(grepl("vfdb", DATABASE)) %>%
  group_by(Treatment, GENE) %>%
  mutate(TPM_log_count = log(TPM + 1)) %>% # Adding 1 to avoid log(0) # Convert the count to relative abundance towards to normal distributions
  select(sample, GENE, TPM_log_count)   # Adding 1 to avoid log(0)


# Compute summary stats
summary_stats <- diversity_filtered %>%
  group_by(Treatment) %>%
  summarise(mean_abun = mean(TPM_log_count),
            se_abun = sd(TPM_log_count) / sqrt(n()))

# Aggregate TPM_log_count values per sample    
geneTable_aggregated <- diversity_filtered %>%
  group_by(sample, Treatment) %>%
  summarise(mean_TPM_log = mean(TPM_log_count, na.rm = TRUE), .groups = "drop") 


# Kruskal–Wallis test (non-parametric alternative to one-way ANOVA)
kruskal_result <- kruskal.test(mean_TPM_log ~ Treatment, data = geneTable_aggregated)
print(kruskal_result)
kruskal_result$p.value
# Passing to ggplot
label_KW <- paste("Kruskal-Wallis p =", format(kruskal_result$p.value, digits = 3))

uniqueGene <- unique(diversity_filtered$GENE)





data_stats = (genetable_normdata %>%
                #filter(DATABASE == "plasmidfinder") %>%
                #filter(DATABASE == "card") %>%
                filter(DATABASE == "vfdb") %>%
                arrange(desc(TPM)) %>%
                mutate(Present = 1) %>%
                pivot_wider(
                  id_cols = c(sample, Treatment),
                  names_from = GENE,
                  values_from = Present,
                  values_fill = list(Present = 0)
                ) %>%
                # ADD RICHNESS CALCULATION HERE
                mutate(Richness = rowSums(pick(3:last_col()))) %>% 
                tbl_summary(
                  by = Treatment,
                  # Include genes (3:13) AND the new Richness column
                  include = c(all_of(3:13), Richness), 
                  type = list(
                    all_of(1:11) ~ "dichotomous",
                    Richness ~ "continuous" # Treat Richness as a number
                  ),
                  value = list(all_of(1:11) ~ 1),
                  statistic = list(
                    all_dichotomous() ~ "{n} / {N} ({p}%)",
                    Richness ~ "{mean} ({sd})" # Show mean richness per treatment
                  ),
                  digits = list(
                    all_dichotomous() ~ c(0, 0, 1),
                    Richness ~ 2
                  ),
                  sort = list(all_of(1) ~ "frequency"),
                  missing = "no"
                )) %>%
  # Kruskal-Wallis for Richness, Chi-sq for genes
  add_p(test = list(Richness ~ "kruskal.test", all_dichotomous() ~ "kruskal.test")) %>%
  add_overall() %>%
  modify_table_body(
    fun = ~ .x %>% 
      dplyr::arrange(
        variable == "Richness",                    # First: Keep Richness at bottom (FALSE < TRUE)
        desc(readr::parse_number(stat_0))          # Second: Sort everything else by count
      )
  ) %>%
  as_gt() %>%
  tab_header(md(glue("**VFs richness by treatment (_{label_KW}_)**")))



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


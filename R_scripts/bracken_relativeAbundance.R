# Purpose: Stacked relative-abundance bar plot (top-50 taxa) and ranked species
#          count bar plot from the Bracken output, faceted by treatment group.
# Input:   data/bracken_arranged.csv  (cols: #perc, tot_all, tot_lvl, taxid, name, samples...)
#          data/chicken_metadata.csv  (cols: sample, Treatment, ...)
#          data/taxa_category.csv     (cols: Taxon, Category)
# Output:  Plots rendered to the active graphics device;
#          interactive HTML saved to output/ (uncomment saveWidget call).

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(purrr)
library(paletteer)
library(plotly)
library(htmlwidgets)
library(stringr)

# Load datasets
bracken_arranged <- read_csv("data/bracken_arranged.csv")
chicken_metadata <- read_csv("data/chicken_metadata.csv")
chicken_metadata <- na.omit(chicken_metadata)
taxa_category    <- read_csv("data/taxa_category.csv")

# Standardise treatment labels
modify_strings <- function(df, column_name, target_string, replacement_string) {
  df %>% mutate(!!sym(column_name) :=
                  str_replace_all(!!sym(column_name), fixed(target_string), replacement_string))
}

chicken_metadata <- modify_strings(chicken_metadata, "Treatment",
                                   "Reference diet (Wheat -soyabean)", "Reference diet")
chicken_metadata <- modify_strings(chicken_metadata, "Treatment",
                                   "Seaweed", "Dulce")

# Pivot longer and filter to target classes
bracken_pivoted <- bracken_arranged %>%
  pivot_longer(cols = 4:21, names_to = "sample", values_to = "count") %>%
  filter(!grepl("root|Homo sapiens|cellular organisms|unclassified|Bacteria|environmental samples", name)) %>%
  filter(grepl("Clostridia|Actinomycetes|Bacilli|Gammaproteobacteria|Alphaproteobacteria|Betaproteobacteria|Bacteroidia|Spirochaetia", name)) %>%
  filter(name != "Clostridiaceae") %>%
  filter(!is.na(name)) %>%
  group_by(sample) %>%
  filter(count > 5) %>%
  mutate(Percentage = (count / sum(count)) * 100) %>%
  ungroup()

# Merge with metadata
bracken_merged <- merge(bracken_pivoted, chicken_metadata, by = "sample")

taxa_category <- taxa_category %>% rename(name = Taxon)
bracken_mergedCategory <- merge(bracken_merged, taxa_category, by = "name")

# Identify top-50 taxa (Unique category only)
top50_taxa <- bracken_mergedCategory %>%
  filter(Category == "Unique") %>%
  filter(name != "B. animalis") %>%
  group_by(name) %>%
  summarise(total = sum(count, na.rm = TRUE)) %>%
  arrange(desc(total)) %>%
  slice_head(n = 50) %>%
  pull(name)

# Collapse non-top-50 taxa into "Others"
bracken_collapsed <- bracken_mergedCategory %>%
  filter(Category == "Unique") %>%
  mutate(name = ifelse(name %in% top50_taxa, name, "Others")) %>%
  group_by(sample) %>%
  mutate(PercTop50 = (count / sum(count)) * 100) %>%
  arrange(sample, desc(PercTop50)) %>%
  ungroup()

bracken_collapsed$name <- factor(
  bracken_collapsed$name,
  levels = bracken_collapsed %>%
    group_by(name) %>%
    summarise(total = sum(PercTop50)) %>%
    arrange(desc(total)) %>%
    pull(name)
)

# Colour palette
taxa       <- as.factor(bracken_merged$name)
t          <- length(levels(taxa))
t_palette  <- paletteer_d("ggsci::default_igv", n = t)

ann_name_col <- bracken_merged %>%
  select(name) %>%
  na.omit() %>%
  distinct(name) %>%
  mutate(color = t_palette)
name_colors <- setNames(as.character(ann_name_col$color), ann_name_col$name)

# Stacked relative-abundance bar plot
ra <- ggplot(bracken_collapsed, aes(x = factor(sample), y = PercTop50, fill = name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = t_palette) +
  labs(x = "Sample Groups by Treatment", y = "Percentage (Relative abundance)",
       fill = "Taxonomic level") +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        text        = element_text(size = 14)) +
  guides(fill = guide_legend(nrow = 10)) +
  facet_wrap(~ Treatment, scales = "free_x", nrow = 1)

ra

# Interactive plotly version
Relative_abundace_plotly <- plotly::ggplotly(ra)
Relative_abundace_plotly

# Uncomment to save as HTML:
# saveWidget(widget = Relative_abundace_plotly,
#            file = "output/relativeAbundance_plot.html",
#            selfcontained = TRUE)

# Ranked species count bar plot (all taxa, no category filter)
df_species <- bracken_merged %>%
  group_by(name) %>%
  summarise(totCount = sum(count), .groups = "drop") %>%
  mutate(PerClass = (totCount / sum(totCount)) * 100) %>%
  filter(totCount > 0) %>%
  arrange(desc(totCount)) %>%
  distinct(totCount, .keep_all = TRUE)

Species_count <- ggplot(df_species, aes(x = reorder(name, totCount), y = totCount, fill = name)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = median(df_species$totCount), linetype = "dashed", color = "black") +
  geom_text(aes(label = paste0(round(totCount / 1000, 1), "K")),
            vjust = 0.5, hjust = -0.05, size = 3.5) +
  scale_fill_manual(values = name_colors) +
  scale_y_log10() +
  labs(title = "", x = "Taxonomy Class Level", y = "Total Count",
       fill = "Taxonomic class") +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        text        = element_text(size = 14)) +
  coord_flip()

Species_count

Totalreads_classLevel <- sum(df_species$totCount)

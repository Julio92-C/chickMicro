# Purpose: Sankey (alluvial) network diagram linking samples -> taxa -> VF genes ->
#          virulence function categories using networkD3.
# Input:   data/abri_kraken2_merged.csv (cols: sample, name, GENE, PRODUCT,
#                                        DATABASE, Treatment, ...)
# Output:  Interactive Sankey diagram in the viewer; uncomment saveWidget to export HTML.

library(plotly)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(networkD3)
library(htmlwidgets)

# Load dataset
abri_kraken2_clean <- read_csv("data/abri_kraken2_merged.csv")

# Filter to VF entries
refactored_df <- abri_kraken2_clean %>%
  filter(DATABASE == "vfdb") %>%
  distinct()

# Extract virulence function category from PRODUCT description
extract_productFunction <- function(string) {
  match  <- regmatches(string, regexpr("\\- ([^\\(]+) \\(", string))
  result <- gsub("\\- | \\(", "", match)
  return(result)
}

refactored_df <- refactored_df %>%
  mutate(Functions = sapply(PRODUCT, extract_productFunction))

# Summary by function
VFs_summary <- refactored_df %>%
  distinct(GENE, .keep_all = TRUE) %>%
  group_by(Functions, GENE) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Functions) %>%
  summarise(Total_count = sum(count, na.rm = TRUE), .groups = "drop")

uniqueVFs  <- length(unique(refactored_df$GENE))
taxaGroup  <- length(unique(refactored_df$name))
uniqueSample <- length(unique(refactored_df$sample))

# Function to add node layer
add_nodes <- function(data, nodes, col_name, color, border) {
  new_nodes <- data %>%
    select(name = {{ col_name }}) %>%
    distinct() %>%
    arrange(name) %>%
    mutate(
      id     = seq(max(nodes$id) + 1, max(nodes$id) + n()),
      color  = color,
      border = border
    )
  bind_rows(nodes, new_nodes)
}

# Build nodes
nodes <- data.frame(
  name = sort(unique(refactored_df$sample)),
  id   = seq(0, length(unique(refactored_df$sample)) - 1)
) %>% mutate(color = "#175709", border = "black")

nodes <- add_nodes(refactored_df, nodes, name,      "#825cdb", "black")
nodes <- add_nodes(refactored_df, nodes, GENE,      "#fc3503", "white")
nodes <- add_nodes(refactored_df, nodes, Functions, "#b5b5b5", "black")

# Function to create links
create_links <- function(data, nodes, source_col, target_col, color) {
  data %>%
    arrange(sample) %>%
    left_join(nodes %>% select(!!source_col := name, source = id), by = source_col) %>%
    left_join(nodes %>% select(!!target_col := name, target = id), by = target_col) %>%
    mutate(color = color) %>%
    select(source, target, color) %>%
    distinct(source, target, .keep_all = TRUE)
}

links0 <- create_links(refactored_df, nodes, "sample",    "name",      "#98ed85")
links1 <- create_links(refactored_df, nodes, "name",      "GENE",      "#4fc1e3")
links2 <- create_links(refactored_df, nodes, "GENE",      "Functions", "#d1d0c24D")

links <- bind_rows(links0, links1, links2) %>%
  distinct(source, target, .keep_all = TRUE) %>%
  arrange(source) %>%
  mutate(value = 1)

# Sankey diagram
sankey_plot <- sankeyNetwork(
  Links      = links,
  Nodes      = nodes,
  Source     = "source",
  Target     = "target",
  Value      = "value",
  NodeID     = "name",
  fontSize   = 14,
  fontFamily = "arial",
  width      = 1500,
  height     = 750,
  sinksRight = TRUE,
  margin     = list(top = 10, right = 10, bottom = 10, left = 10)
)

sankey_plot

# Uncomment to save as HTML:
# saveWidget(sankey_plot, "output/VFsProfile_sankeyNetwork.html")

# Load library
library(tibble)
library(dplyr)
library(igraph)
library(ggraph)
library(readr)
library(stringr)
library(vegan)
library(tidygraph) # Assuming your graph 'g' is a tbl_graph object or similar
library(ggrepel)   # For geom_node_text(repel = TRUE)
library(ggforce)   # For geom_mark_hull (if you want to mark clusters)
library(funrar) # transform an abundance matrix into a relative abundance matrix


# Load dataset
microbiome_df <- read_csv("../data/abri_kraken2Bracken_merged.csv")



# # Filter by Species and counts greater than 30
df_species <- microbiome_df %>%
  group_by(name) %>%
  filter(!grepl("Terrabacteria group|Bacteroidota/Chlorobiota group|FCB group", name)) %>%
  filter(str_count(name, "\\S+") > 1) %>%
  mutate(
    # 1. Abbreviate Genus
    name = str_replace(name, "^(\\w)\\w+\\s(\\w+)", "\\1. \\2"),
    
    # 2. Truncate if too long
    name = if_else(
      str_length(name) > 15,
      str_replace(name, "^((?:\\S+\\s+){2}).*", "\\1"),
      name
    ), 
    
    # 3. Clean up trailing white space
    name = str_trim(name, side = "right") 
  ) %>%
  mutate(name = case_when(
    grepl("\\[Clostridium\\] asparagiforme ", name) ~ "C. asparagiforme",
    grepl("\\[Clostridium\\] scindens", name) ~ "C. scindens",
    grepl("\\[Ruminococcus\\] lactaris", name) ~ "R. lactaris ",
    grepl("\\[Ruminococcus\\] torques", name) ~ "R. torques",
    grepl("\\[Clostridium\\] innocuum", name) ~ "C. innocuum",
    TRUE ~ name
  ),
  name = ifelse(grepl("B. pullorum subsp. gallinarum", name), "B. pullorum", name)) %>%
  mutate(GENE = gsub("Bifidobacterium_adolescentis_rpoB_mutants_conferring_resistance_to_rifampicin", "rpoB_mutants", GENE)) %>%
  mutate(GENE = gsub("Bifidobacterium_bifidum_ileS_conferring_resistance_to_mupirocin", "ileS", GENE)) %>%
  #summarise(count = sum(count), .groups = "drop") %>%
  filter(sampleCount > 10) %>%
  distinct() %>%
  rename(taxa = name)



# Unique species
uniqueTaxa <- unique(df_species$taxa)

# Define taxa groups
lab_taxa <- c("L. salivarius", "L. phocaeense", "L. crispatus", "L. asacharolyticus", "B. pseudolongum", "B. pseudolongum PV8-2", "B. pullorum")


scfa_taxa <- c("A. finegoldii", "C. scindens", "A. finegoldii DSM 17242", "F. prausnitzii", "F. plautii", "B. obeum", "B. obeum ATCC 29174", "B. obeum A2-162", "B. argi", "B. wexlerae", "B. wexlerae DSM 19850", "R. lactaris", "A. hadrus", "B. hansenii", "B. hansenii DSM 20583", "M. formatexigens", "M. formatexigens DSM 14469", "B. pseudococcoides", "F. sp. I3-3-89", "B. parvula", "A. hallii", "B. producta", "B. animalis", "B. animalis", "D. longicatena", "R. intestinalis", "B. longum")


pathon_taxa <- c("C. sp. M62/1", "C. sp. C1", "C. bilis", "C. difficile", "C. asparagiforme", "E. coli", "E. bolteae", "E. asparagiformis", "E. callanderi", "E. incertae sedis", "S. intestinalis", "E. faecium", "C. perfringens")

mucgly_taxa <- c("R. torques", "M. gnavus", "u. Subdoligranulum", "S. variabile", "C. sp. M62/1")


hydrog_taxa <- c("B. hydrogenotrophica", "B. hydrogenotrophica DSM 10507")


# Add functional group column
df_species <- df_species %>%
  mutate(Function_group = case_when(
    taxa %in% lab_taxa ~ "Lactic acid bacteria/probiotics",
    taxa %in% scfa_taxa ~ "Primary fermenters / SCFA producers",
    taxa %in% mucgly_taxa ~ "Mucin/glycan degraders and niche specialists",
    taxa %in% hydrog_taxa ~ "Hydrogenotrophs / cross‑feeders",
    taxa %in% pathon_taxa ~ "Proteolytic/opportunistic taxa and pathobionts",
    TRUE ~ "Other/unknown"
  )) 


# 1. Prepare Edges (including new links between Samples and Treatments)

# Weighted Sample → Taxa edges
sample_taxa_edges <- df_species %>%
  select(from = sample, to = taxa, weight = sampleCount)

# Unweighted Taxa → Gene edges (can be binary or frequency-based later)
taxa_gene_edges <- df_species %>%
  select(from = taxa, to = GENE) %>%
  mutate(weight = 1)

# Combine edges
edges <- bind_rows(sample_taxa_edges, taxa_gene_edges)

# 2. Prepare Node Metadata (including Treatment Groups as their own nodes)

# Sample metadata
sample_meta <- df_species %>%
  select(sample, Treatment) %>%
  distinct()

# Taxa metadata
taxa_meta <- df_species %>%
  select(taxa, Function_group) %>%
  distinct()

# Gene metadata
gene_meta <- df_species %>%
  group_by(GENE) %>%
  select(GENE, DATABASE) %>%
  rename(gene = GENE) %>%
  rename(geneCategory = DATABASE) %>%
  distinct()  %>% 
  mutate(geneCategory = case_when(
    grepl("card", geneCategory) ~ "AMR",
    grepl("vfdb", geneCategory) ~ "VFs",
    grepl("plasmidfinder", geneCategory) ~ "MGEs",
    TRUE ~ geneCategory
  ))
  

# Build node list
nodes <- unique(c(edges$from, edges$to)) %>%
  tibble(name = .) %>%
  mutate(type = case_when(
    name %in% df_species$sample ~ "sample",
    name %in% df_species$taxa   ~ "taxa",
    name %in% df_species$GENE   ~ "gene"
  )) %>%
  left_join(sample_meta, by = c("name" = "sample")) %>%
  left_join(taxa_meta, by = c("name" = "taxa")) %>%
  left_join(gene_meta,   by = c("name" = "gene")) %>%
  distinct(name, .keep_all = TRUE) # remove duplicates

# 3. (Optional) Your Bray-Curtis Clustering analysis
# This analysis is valid for grouping SAMPLES, but is separate from the network itself
# and just adds metadata to the node table.

# Create a sample × taxa matrix from your microbiome_df
abundance_matrix <- df_species %>%
  select(sample, taxa, sampleCount) %>%
  # mutate(
  #   count = (count / sum(count)) * 100 # Convert the count to relative abundance towards to normal distributions
  # ) %>%
  tidyr::pivot_wider(
    names_from = taxa,
    values_from = sampleCount,
    values_fill = 0,
    values_fn = sum # Use `sum` to aggregate duplicate values
  ) %>%
  # distinct(sample, .keep_all = TRUE) %>%
  column_to_rownames("sample")

# Convert the data frame to a matrix
abundance_matrix <- as.matrix(abundance_matrix)

# transform an abundance matrix into a relative abundance matrix
relative_abundance_matrix <- abundance_matrix %>%
  make_relative()

# Compute Bray-Curtis distance
dist_matrix <- vegdist(relative_abundance_matrix, method = "bray")

# Hierarchical clustering
hc <- hclust(dist_matrix, method = "average")

# Cut tree into clusters (e.g. 2 groups)
sample_clusters <- cutree(hc, k = 4)

# Add cluster info to node table
nodes$distance_cluster <- sample_clusters[match(nodes$name, names(sample_clusters))]


## 4. Prepare tables for Gephi Export

# Clean edge table for Gephi
edges_gephi <- edges %>%
  rename(Source = from, Target = to, Weight = weight)

summary(edges_gephi)

# Clean node table for Gephi
nodes_gephi <- nodes %>%
  rename(Id = name) %>%
  select(Id, type, Treatment, geneCategory, Function_group, distance_cluster)


# Replace NA values with an empty string ""
# Replace NA values with an empty string ""
nodes_gephi$Treatment[is.na(nodes_gephi$Treatment)] <- ""
nodes_gephi$geneCategory[is.na(nodes_gephi$geneCategory)] <- ""
nodes_gephi$Function_group[is.na(nodes_gephi$Function_group)] <- ""
nodes_gephi$distance_cluster[is.na(nodes_gephi$distance_cluster)] <- ""


# Class 
summary(nodes_gephi)

# Export node and edges to CSV files
write.csv(edges_gephi, "../data/gephi_edges3.csv", row.names = FALSE)
write.csv(nodes_gephi, "../data/gephi_nodes3.csv", row.names = FALSE)


# Create graph object
g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

# --- Start of your modified plotting code ---

# Plot with gene nodes colored by GeneCategory
ggraph(g, layout = "fr",
       niter = 1000, # Increase iterations for more stable layout, potentially better spread
       area = vcount(g)^2 * 1.5) + # Increase 'area' parameter for Fruchterman-Reingold
  # to potentially spread nodes further apart
  
  geom_edge_link(aes(width = weight), 
                 alpha = 0.2) + # Reduce edge opacity to thin out dense areas
  # (original was 0.4, try 0.2 or even lower like 0.1)
  
  geom_node_point(aes(
    shape = type,
    color = case_when(
      type == "sample" ~ as.factor(distance_cluster),
      type == "gene"   ~ geneCategory,
      TRUE             ~ type
    )
  ), 
  size = 5,
  alpha = 0.7) + # Add alpha to node points for transparency, especially for overlap
  # Adjust 'alpha' to your preference (e.g., 0.7-0.9)
  
  geom_node_text(aes(label = name), 
                 repel = TRUE, 
                 size = 3, 
                 max.overlaps = Inf,
                 segment.color = 'grey50', # Add a segment line for clarity when labels are repelled
                 min.segment.length = 0.5) + # Only draw segments if they are long enough
  
  scale_edge_width(range = c(0.5, 3)) +
  scale_color_manual(values = c(
    "1" = "#1f77b4", "2" = "#ff7f0e", "3" = "yellow", "4"="green", # Sample clusters
    AMR         = "#d62728",
    VFs         = "#9467bd",
    MGEs        = "#2ca02c",
    taxa        = "#8c564b"
  )) +
  scale_shape_manual(values = c(sample = 16, taxa = 17, gene = 15)) +
  
  # --- Optional: Add geom_mark_hull to highlight clusters if needed ---
  # Assuming 'distance_cluster' or 'geneCategory' defines your central clusters
  # You might want to filter this to specific clusters if the entire center is a mix
  # For example, to highlight 'gene' nodes by their 'geneCategory':
  # geom_mark_hull(aes(x = x, y = y, fill = geneCategory, filter = type == "gene"),
  #                alpha = 0.1, show.legend = FALSE) +
  
  theme_void() +
  theme(legend.position = "none")

# --- End of your modified plotting code ---


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

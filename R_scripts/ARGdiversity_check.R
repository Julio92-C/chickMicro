# Purpose: Exploratory diversity checks for ARG profiles: rank-abundance curves,
#          top-N dominance, prevalence, richness sensitivity, Shannon diversity,
#          and a heatmap of the top 20 ARGs.
# Input:   data/genetable_normdata.csv (cols: sample, GENE, DATABASE, TPM, Treatment, ...)
# Output:  Plots and summary output rendered to the active graphics device / console.

library(tidyverse)
library(vegan)

# Load dataset
genetable_normdata <- read_csv("data/genetable_normdata.csv")

# Metadata
meta <- genetable_normdata %>%
  distinct(sample, .keep_all = TRUE) %>%
  select(sample, Treatment) %>%
  column_to_rownames("sample")

# Build log-transformed ARG abundance matrix (GENE x sample)
geneTable_filtered <- genetable_normdata %>%
  filter(grepl("card", DATABASE)) %>%
  mutate(TPM_log_count = log(TPM + 1)) %>%
  select(sample, GENE, TPM_log_count) %>%
  tidyr::pivot_wider(
    names_from  = sample,
    values_from = TPM_log_count,
    values_fill = 0,
    values_fn   = sum
  ) %>%
  distinct(GENE, .keep_all = TRUE) %>%
  column_to_rownames("GENE")

tpm_mat <- as.matrix(geneTable_filtered)

# 1. Rank-abundance curve for the first sample
sample_id <- colnames(tpm_mat)[1]
ranked    <- sort(tpm_mat[, sample_id], decreasing = TRUE)
cumfrac   <- cumsum(ranked) / sum(ranked)
plot(cumfrac, type = "l", main = sample_id, xlab = "Ranked ARG", ylab = "Cumulative fraction")

# 2. Top-N cumulative fraction across treatment groups
topN         <- 5
top_frac_by_sample <- apply(tpm_mat, 2, function(x) sum(sort(x, decreasing = TRUE)[1:topN]) / sum(x))
df_top <- tibble(sample = names(top_frac_by_sample), top_frac = top_frac_by_sample) %>%
  left_join(meta %>% rownames_to_column("sample"), by = "sample")
ggplot(df_top, aes(x = Treatment, y = top_frac)) + geom_boxplot() + geom_jitter(width = 0.2)

# 3. Per-feature mean TPM and prevalence
feat_stats <- tibble(
  ARG        = rownames(tpm_mat),
  mean_TPM   = rowMeans(tpm_mat),
  prevalence = rowSums(tpm_mat > 0)
)
head(arrange(feat_stats, desc(mean_TPM)), 20)

# 4. Richness sensitivity to detection threshold
thresholds <- c(0, 0.1, 1, 5)
richness_by_thresh <- map_dfr(thresholds, function(th) {
  r <- colSums(tpm_mat > th)
  tibble(threshold = th, sample = names(r), richness = r)
})
richness_by_thresh <- richness_by_thresh %>%
  left_join(meta %>% rownames_to_column("sample"), by = "sample")
ggplot(richness_by_thresh, aes(x = factor(threshold), y = richness, fill = Treatment)) +
  geom_boxplot()

# 5. Shannon diversity
shannon <- apply(tpm_mat, 2, function(x) {
  p <- x / sum(x)
  diversity(p, index = "shannon")
})
richness <- colSums(tpm_mat > 0)
alpha_df <- tibble(
  sample   = colnames(tpm_mat),
  Treatment = meta[colnames(tpm_mat), "Treatment"],
  richness  = richness,
  shannon   = shannon
)

# 6. Heatmap of top 20 ARGs by mean TPM
top20 <- feat_stats %>%
  arrange(desc(mean_TPM)) %>%
  slice_head(n = 20) %>%
  pull(ARG)

library(pheatmap)
pheatmap(log10(tpm_mat[top20, ] + 1), annotation_col = meta, cluster_rows = TRUE, cluster_cols = FALSE)

# 7. Outlier detection: samples with extreme total TPM
total_ARG_TPM <- colSums(tpm_mat)
boxplot(total_ARG_TPM ~ meta$Treatment)
which(total_ARG_TPM > quantile(total_ARG_TPM, 0.95))

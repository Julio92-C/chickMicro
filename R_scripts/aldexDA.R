# Purpose: Basic ALDEx2 differential abundance workflow: overall Kruskal-Wallis
#          (all three treatments) plus pairwise t-tests between each pair of
#          treatment groups, with Bland-Altman, effect, and volcano plots.
# Input:   data/Bracken_SpeciesCounts.csv (cols: name, sample columns...)
#          data/chicken_metadata1.csv     (cols: sample, Treatment, ...)
# Output:  ALDEx2 diagnostic plots rendered to the active graphics device;
#          results written to output/aldex2_*_results.tsv (uncomment as needed).

library(ALDEx2)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(stringr)

# Load inputs
Bracken_SpeciesCounts <- read_csv("data/Bracken_SpeciesCounts.csv")
chicken_metadata      <- read_csv("data/chicken_metadata1.csv")

# Move names to row names
counts <- Bracken_SpeciesCounts %>%
  distinct(name, .keep_all = TRUE) %>%
  column_to_rownames("name")

# Move sample to row names and arrange by treatment
meta <- chicken_metadata %>%
  arrange(Treatment) %>%
  column_to_rownames("sample")

# Sanity checks
counts <- counts[, rownames(meta)]   # reorder columns to match metadata
stopifnot(identical(colnames(counts), rownames(meta)))

# Filter low-prevalence taxa (keep taxa present in >= 2 samples)
keep     <- rowSums(counts > 0) >= 2
counts_f <- counts[keep, ]

# Run aldex DAA across all treatments
results <- aldex(reads = counts_f, conditions = meta$Treatment, mc.samples = 128,
                 test = "kw", effect = FALSE, CI = FALSE,
                 include.sample.summary = FALSE,
                 verbose = FALSE, paired.test = FALSE,
                 denom = "all", iterate = FALSE, gamma = NULL)

# Pairwise comparison — Dulce vs Reference diet
counts_f_subset <- counts_f[, 1:12]
meta_subset     <- meta %>% slice(1:12)

# Pairwise comparison — Reference diet vs Soyabean meal
counts_f_subset <- counts_f[, 7:18]
meta_subset     <- meta %>% slice(7:18)

# Pairwise comparison — Dulce vs Soyabean meal
counts_f_subset <- counts_f[, c(1:6, 13:18)]
meta_subset     <- meta %>% slice(c(1:6, 13:18))

# Run aldex DAA (pairwise)
results <- aldex(reads = counts_f_subset, conditions = meta_subset$Treatment,
                 mc.samples = 128,
                 test = "t", effect = TRUE, CI = FALSE,
                 include.sample.summary = FALSE,
                 verbose = FALSE, paired.test = FALSE,
                 denom = "all", iterate = FALSE, gamma = NULL)

# Diagnostic plots
par(mfrow = c(1, 3))
aldex.plot(results, type = "MA", test = "welch",
           xlab = "Log-ratio abundance", ylab = "Difference",
           main = "Bland-Altman plot")
aldex.plot(results, type = "MW", test = "welch",
           xlab = "Dispersion", ylab = "Difference",
           main = "Effect plot")
aldex.plot(results, type = "volcano", test = "welch",
           xlab = "Difference", ylab = "-1(log10(q))",
           main = "Volcano plot")

# Top results by effect size
top_by_effect <- results[order(-abs(results$effect)), ]
head(top_by_effect, 10)

# Filter candidates (example thresholds)
candidates <- subset(results, abs(effect) >= 1 & overlap <= 0.10)
candidates[order(-abs(candidates$effect)), ]

# Boxplot of raw counts for a taxon
taxon <- "B. hydrogenotrophica"
boxplot(as.numeric(counts[taxon, ]) ~ meta$Treatment,
        ylab = "Raw counts", xlab = "Treatment", main = taxon)
stripchart(as.numeric(counts[taxon, ]) ~ meta$Treatment,
           add = TRUE, vertical = TRUE, pch = 16, col = "blue")

# Volcano-style plot (effect vs -log10 q)
plot(results$effect, -log10(results$wi.eBH), pch = 20,
     xlab = "ALDEx2 effect (CLR diff)", ylab = "-log10(wilcoxon BH p)")
abline(v = c(-1, 1), col = "grey", lty = 2)
abline(h = -log10(0.05), col = "red", lty = 2)

# Write out results — uncomment the relevant comparison before running:
# write.table(results, file = "output/aldex2_results.tsv",
#             sep = "\t", quote = FALSE, col.names = NA)
# write.table(results, file = "output/Dulce_vs_Referencediet_aldex2_results.tsv",
#             sep = "\t", quote = FALSE, col.names = NA)
# write.table(results, file = "output/ReferenceDiet_vs_SoyabeanMeal_aldex2_results.tsv",
#             sep = "\t", quote = FALSE, col.names = NA)
# write.table(results, file = "output/Dulce_vs_SoyabeanMeal_aldex2_results.tsv",
#             sep = "\t", quote = FALSE, col.names = NA)

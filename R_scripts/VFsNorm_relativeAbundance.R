# Purpose: Stacked relative-abundance bar plots and ranked bar charts for VF
#          functional categories from TPM-normalised data, with Kruskal-Wallis tests.
# Input:   data/genetable_normdata.csv (cols: sample, GENE, PRODUCT, DATABASE,
#                                       TPM, Treatment, ...)
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

# Extract virulence function category from PRODUCT description
extract_productFunction <- function(string) {
  match  <- regmatches(string, regexpr("\\- ([^\\(]+) \\(", string))
  result <- gsub("\\- | \\(", "", match)
  return(result)
}

# Compute per-sample percentage of VF reads
genetable_filtered <- genetable_normdata %>%
  filter(DATABASE == "vfdb") %>%
  mutate(Functions = sapply(PRODUCT, extract_productFunction)) %>%
  select(sample, GENE, PRODUCT, Functions, TPM, Treatment) %>%
  group_by(sample) %>%
  mutate(Percentage = (TPM / sum(TPM)) * 100) %>%
  ungroup()

uniqueVFs      <- length(unique(genetable_filtered$GENE))
uniqueFunction <- length(unique(genetable_filtered$Functions))

Functions_f <- as.factor(genetable_filtered$Functions)
f           <- length(levels(Functions_f))
f_palette   <- paletteer_d("ggsci::default_nejm", n = f)

# Stacked relative-abundance bar plot
ra <- ggplot(genetable_filtered, aes(x = factor(sample), y = Percentage, fill = Functions)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = f_palette) +
  labs(x = "Sample Groups by Treatment", y = "Percentage (Relative abundance)") +
  theme_classic() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14)) +
  guides(fill = guide_legend(nrow = 10)) +
  labs(fill = "Virulent functions") +
  facet_wrap(~ Treatment, scales = "free_x", nrow = 1)

ra

# Ranked total TPM per virulence function
vfdb_totals <- genetable_filtered %>%
  group_by(Functions) %>%
  summarise(Total_TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total_TPM))

vfdb_totals[8, 1] <- "Antimicrobial activity/CA*"

VFs_count <- ggplot(vfdb_totals, aes(x = reorder(Functions, Total_TPM),
                                      y = Total_TPM, fill = Functions)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = median(vfdb_totals$Total_TPM), linetype = "dashed", color = "black") +
  geom_text(aes(label = paste0(round(Total_TPM / 1000, 1), "K")),
            vjust = 0.5, hjust = 1.0, size = 3.5) +
  scale_fill_manual(values = f_palette) +
  scale_y_log10() +
  labs(title = "", x = "Virulent Functions", y = "Total Count") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        text = element_text(size = 14)) +
  coord_flip()

VFs_count

# Subset to one functional category and run Kruskal-Wallis + violin plot
vfdb_filtered <- genetable_filtered %>%
  filter(Functions == "Nutritional/Metabolic factor") %>%
  mutate(log_TPM = log(TPM + 1)) %>%
  arrange(desc(log_TPM))

summary_stats <- vfdb_filtered %>%
  group_by(Treatment) %>%
  summarise(mean_count = mean(log_TPM),
            se_count   = sd(log_TPM) / sqrt(n()))

classTable_aggregated <- vfdb_filtered %>%
  group_by(sample, Treatment) %>%
  summarise(mean_log_TPM = mean(log_TPM, na.rm = TRUE), .groups = "drop")

kruskal_result <- kruskal.test(mean_log_TPM ~ Treatment, data = classTable_aggregated)
print(kruskal_result)
label_KW <- paste("Kruskal-Wallis p =", format(kruskal_result$p.value, digits = 2))

my_comparisons <- list(
  c("Reference diet", "Dulce"),
  c("Soyabean meal",  "Reference diet"),
  c("Dulce",          "Soyabean meal")
)

ggplot(vfdb_filtered, aes(x = Treatment, y = log_TPM, fill = Treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                               "Dulce"          = "#2ca02c",
                               "Soyabean meal"  = "#ff7f0e")) +
  labs(x = "Sample groups", y = "VFs abundance (log10(TPM))") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(max(vfdb_filtered$log_TPM) + 0.5,
                                 max(vfdb_filtered$log_TPM) + 1.3,
                                 max(vfdb_filtered$log_TPM) + 1.8),
                     method = "wilcox.test",
                     label  = "p.signif") +
  annotate("text", x = 2, y = max(vfdb_filtered$log_TPM) + 3,
           label = label_KW, size = 4, color = "black")

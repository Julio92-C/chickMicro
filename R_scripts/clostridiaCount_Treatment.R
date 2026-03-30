# Load libraries
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



# Import Bracken arrange report 
bracken_arranged <- read_csv("../data/bracken_arranged.csv")
View(bracken_arranged)
colnames(bracken_arranged)

# Load metadata
chicken_metadata <- read_csv("../data/chicken_metadata.csv")
chicken_metadata <- na.omit(chicken_metadata)
View(chicken_metadata)  


# Pivot longer Bracken report
bracken_pivoted <- bracken_arranged %>%
  pivot_longer(cols = 4:21, names_to = "sample", values_to = "count") %>%
  filter(!grepl("root|Homo sapiens|cellular organisms|unclassified|Bacteria|environmental samples", name)) %>%
  filter(!is.na(name)) %>%
  group_by(sample) %>%
  filter(count > 1) %>%
  mutate(Percentage = (count / sum(count)) * 100) %>%
  ungroup()

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
chicken_metadata <- modify_strings(chicken_metadata, "Treatment", 
                                   "Seaweed", 
                                   "Dulce")

# Merge bracken_arranged and metadata dataframe
bracken_merged <- merge(bracken_pivoted, chicken_metadata, by ="sample", all.x = TRUE)


# Filter the dataset by Class
class_data <- bracken_merged %>%
  #filter(name == "Clostridia")
  #filter(name == "Actinomycetes")
  #filter(name == "Bacilli")
  #filter(name == "Gammaproteobacteria")
  #filter(name== "Alphaproteobacteria")
  #filter(name == "Betaproteobacteria")
  #filter(name == "Bacteroidota")
  #filter(name == "Spirochaetia")
  #filter(name == "Alistipes finegoldii") %>%
  #filter(name == "[Clostridium] scindens") %>%
  filter(name == "Enterobacteriaceae") %>%
  # filter(name == "Ligilactobacillus salivarius") %>%
  # filter(name == "Clostridium perfringens") %>%
  #filter(name == "Clostridioides difficile") %>%
  #filter(name == "Enterococcus faecium") %>%
  #filter(name == "Escherichia coli") %>%
  mutate(log_count = log(count + 1))


# Compute summary stats
summary_stats <- class_data %>%
  group_by(Treatment) %>%
  summarise(mean_count = mean(count),
            se_count = sd(count) / sqrt(n()))

# Aggregate log_count values per sample    
classTable_aggregated <- class_data %>%
  group_by(sample, Treatment) %>%
  summarise(mean_log_count = mean(log_count, na.rm = TRUE), .groups = "drop")



# Perform ANOVA
# anova_result <- aov(mean_log_count ~ Treatment, data = classTable_aggregated)
# anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
# print(anova_p)



# Kruskal–Wallis test (non-parametric alternative to one-way ANOVA)
kruskal_result <- kruskal.test(mean_log_count ~ Treatment, data = classTable_aggregated)
print(kruskal_result)
kruskal_result$p.value
# Passing to ggplot
label_KW <- paste("Kruskal-Wallis p =", format(kruskal_result$p.value, digits = 2))



# Compute mean and SE per sample
sample_stats <- class_data %>%
  group_by(sample) %>%
  summarise(mean_count = mean(count),
            se_count = sd(count) / sqrt(n()),
            Treatment = first(Treatment)) # retain Treatment info


# Class count
Species_count <- ggplot(class_data, aes(x = sample, y = count, fill = Treatment)) +
  geom_bar(stat = "identity") +  # Use stat = "identity" for pre-calculated counts
  geom_text(aes(label = round(count/1,1)), vjust = -0.6, size = 3) + # Explicitly define y for labels
  geom_hline(aes(yintercept = mean(count)), linetype = "dashed", color = "red") +
  annotate("text", x = Inf, y = mean(class_data$count), label = "", hjust = 1.1, vjust = -0.5, color = "red", size = 3)+
  geom_line(aes(group = 1), color = "black", linewidth = 1) +
  geom_errorbar(data = sample_stats,
                aes(x = sample, ymin = mean_count - se_count, ymax = mean_count + se_count),
                width = 0.2, color = "white", inherit.aes = FALSE) +
  geom_point(data = sample_stats, aes(x = sample, y = mean_count), color = "red", size = 2, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                               "Dulce" = "#2ca02c",
                               "Soyabean meal" = "#ff7f0e")) +
  labs(
    title = bquote(italic("Enterobacteriaceae") ~ "Count per Treatment (Kruskal-Wallis p =" ~ .(signif(kruskal_result$p.value, 3)) ~ ")"),
    x = "Sample groups",
    y = "Count"
  ) +
  theme_classic() +  # Apply classic theme
  theme(legend.position = "none", text = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust = 0.5)
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        text = element_text(size = 14)) + # Set text size
        facet_wrap(~ Treatment, scales = "free_x", nrow = 1) 

# Plot species count
Species_count



# Define all 3 pairwise comparisons
my_comparisons <- list( 
  c("Reference diet", "Dulce"),
  c("Soyabean meal", "Reference diet"),
  c("Dulce", "Soyabean meal")
  
)


# Generate the violin plot with statistical significance
ggplot(class_data, aes(x = Treatment, y = log_count, fill = Treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  # Jittered points for individual samples
  #geom_jitter(width = 0.1, size = 2, alpha = 0.3) +
  scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                               "Dulce" = "#2ca02c",
                               "Soyabean meal" = "#ff7f0e")) +
  labs(x = "Sample groups", y = "Enterobacteriaceae (log10(Cont))") +
  theme_classic() +
  theme(legend.position = "none", text = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust = 0.5)
  ) +
  # Add the pairwise comparisons using Wilcoxon test
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(
                       max(class_data$log_count) + 0.5,  # Control_4 vs Control_5
                       max(class_data$log_count) + 1.3,  # Control_4 vs Seaweed_4
                       max(class_data$log_count) + 1.8  # Control_5 vs Seaweed_5
                     ),
                     method = "wilcox.test",
                     label = "p.signif") + # Shows asterisks (*, **, ***)
  #geom_text(aes(label = Letters), vjust = -0.5, size = 5) +
  annotate("text", x = 2, y = max(class_data$log_count) + 3,
           label = label_KW,
           size = 4, color = "black") 
coord_flip()


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
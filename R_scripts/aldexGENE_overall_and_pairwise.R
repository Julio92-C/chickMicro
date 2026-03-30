#!/usr/bin/env Rscript
# Purpose: ALDEx2 overall (Kruskal-Wallis) + pairwise (t-test) automated pipeline
#          for gene-level differential abundance analysis, with publication-ready
#          diagnostic plots (MA, MW, volcano, bar, violin) written to PDF and CSV.
# Input:   data/geneTable_rawCounts.csv (cols: name, sample columns...)
#          data/chicken_metadata1.csv   (cols: sample, Treatment, ...)
# Output:  output/aldex_output_*.tsv and *.csv; PDF plot files in working directory.

library(ALDEx2)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)
library(ggpubr)
library(plotly)

# -------------------------
# User inputs / parameters
# -------------------------
counts_path <- "data/geneTable_rawCounts.csv"
meta_path   <- "data/chicken_metadata1.csv"

out_prefix          <- "output/aldex_output"
mc.samples_overall  <- 128
mc.samples_pairwise <- 128
denom_choice        <- "all"   # try "iqlr" for sensitivity
set.seed(12345)

# -------------------------
# Read inputs
# -------------------------
message("Reading input files...")
geneTable_rawCounts <- read_csv(counts_path, show_col_types = FALSE)
chicken_metadata    <- read_csv(meta_path,   show_col_types = FALSE)

# Ensure expected column names exist
if (!"name"      %in% colnames(geneTable_rawCounts)) stop("Counts CSV must contain a 'name' column with gene.")
if (!"sample"    %in% colnames(chicken_metadata))    stop("Metadata CSV must contain a 'sample' column.")
if (!"Treatment" %in% colnames(chicken_metadata))    stop("Metadata CSV must contain a 'Treatment' column.")

# Move names to rownames (gene x samples)
counts <- geneTable_rawCounts %>%
  distinct(name, .keep_all = TRUE) %>%
  column_to_rownames("name")

# Move sample to rownames and arrange by Treatment
meta <- chicken_metadata %>%
  arrange(Treatment) %>%
  column_to_rownames("sample")

# -------------------------
# Sanity checks and alignment
# -------------------------
if (!all(rownames(meta) %in% colnames(counts))) {
  if (all(rownames(meta) %in% rownames(counts))) {
    message("Counts appear to be samples x gene; transposing counts.")
    counts <- t(counts)
  } else {
    stop("Unable to align samples between counts and metadata. Check sample IDs in both files.")
  }
}

counts <- counts[, rownames(meta), drop = FALSE]
stopifnot(identical(colnames(counts), rownames(meta)))
message("Counts and metadata aligned: ", nrow(counts), " gene x ", ncol(counts), " samples.")

# Filter low-prevalence genes
min_samples_present <- 2
keep     <- rowSums(as.matrix(counts) > 0) >= min_samples_present
counts_f <- counts[keep, , drop = FALSE]
message("Filtered gene: kept ", nrow(counts_f), " gene (present in >= ", min_samples_present, " samples).")

# -------------------------
# Overall (omnibus) ALDEx2: Kruskal-Wallis
# -------------------------
message("Running overall (Kruskal-Wallis) ALDEx2 test...")
overall_res <- aldex(reads      = counts_f,
                     conditions = meta$Treatment,
                     mc.samples = mc.samples_overall,
                     test       = "kw",
                     effect     = FALSE,
                     CI         = FALSE,
                     include.sample.summary = FALSE,
                     verbose    = FALSE,
                     paired.test = FALSE,
                     denom      = denom_choice,
                     iterate    = FALSE,
                     gamma      = NULL)

overall_outfile <- paste0(out_prefix, "_overall_kw.tsv")
write.table(as.data.frame(overall_res), file = overall_outfile,
            sep = "\t", quote = FALSE, col.names = NA)
message("Overall results written to: ", overall_outfile)

# -------------------------
# Pairwise comparisons (label-based)
# -------------------------
comparisons <- list(
  Dulce_vs_Reference   = c("Dulce", "Reference diet"),
  Reference_vs_Soybean = c("Reference diet", "Soyabean meal"),
  Dulce_vs_Soybean     = c("Dulce", "Soyabean meal")
)

aldex_args <- list(mc.samples = mc.samples_pairwise,
                   test       = "t",
                   effect     = TRUE,
                   CI         = FALSE,
                   include.sample.summary = FALSE,
                   verbose    = FALSE,
                   paired.test = FALSE,
                   denom      = denom_choice,
                   iterate    = FALSE,
                   gamma      = NULL)

run_aldex_on_subset <- function(counts_all, meta_all, groups_vec, comp_name, aldex_args) {
  keep_samples <- rownames(meta_all)[meta_all$Treatment %in% groups_vec]
  if (length(keep_samples) < 2) stop("Not enough samples for comparison: ", comp_name)
  meta_sub   <- meta_all[keep_samples, , drop = FALSE]
  counts_sub <- counts_all[, keep_samples, drop = FALSE]
  counts_sub <- counts_sub[, rownames(meta_sub), drop = FALSE]
  stopifnot(identical(colnames(counts_sub), rownames(meta_sub)))
  message("  ", comp_name, ": running on ", ncol(counts_sub), " samples (",
          paste(groups_vec, collapse = " vs "), ")")
  res    <- do.call(aldex, c(list(reads = counts_sub, conditions = meta_sub$Treatment), aldex_args))
  res_df <- as.data.frame(res)
  res_df$GENE <- rownames(res_df)
  prefixed <- res_df %>%
    relocate(GENE) %>%
    rename_with(.cols = -GENE, .fn = ~ paste0(comp_name, "_", .))
  return(prefixed)
}

results_list <- list()
for (comp in names(comparisons)) {
  groups_vec           <- comparisons[[comp]]
  results_list[[comp]] <- run_aldex_on_subset(counts_f, meta, groups_vec, comp, aldex_args)
  write.table(results_list[[comp]],
              file = paste0(out_prefix, "_", comp, ".tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved pairwise results: ", paste0(out_prefix, "_", comp, ".tsv"))
}

# -------------------------
# Merge pairwise results
# -------------------------
message("Merging pairwise results...")
merged_pairwise <- reduce(results_list, full_join, by = "GENE")

overall_df            <- as.data.frame(overall_res)
overall_df$GENE       <- rownames(overall_df)
keep_overall_cols     <- intersect(c("kw.ep", "kw.eBH", "kw.effect"), colnames(overall_df))
if (length(keep_overall_cols) == 0) keep_overall_cols <- setdiff(colnames(overall_df), "GENE")
overall_subset <- overall_df %>% select(GENE, all_of(keep_overall_cols))

merged_all <- left_join(merged_pairwise, overall_subset, by = "GENE")

merged_outfile <- paste0(out_prefix, "_merged_pairwise_with_overall.tsv")
write.table(merged_all, file = merged_outfile,
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Merged pairwise + overall results written to: ", merged_outfile)

# Compute mean and SE per GENE x Treatment from raw counts
rawCount_stats <- geneTable_rawCounts %>%
  distinct(name, .keep_all = TRUE) %>%
  pivot_longer(cols = -name, names_to = "sample", values_to = "count") %>%
  inner_join(chicken_metadata, by = "sample") %>%
  group_by(name, Treatment) %>%
  summarise(
    mean_count = mean(count, na.rm = TRUE),
    se_count   = ifelse(n() > 1, sd(count, na.rm = TRUE) / sqrt(n()), NA_real_),
    n_samples  = n(),
    .groups    = "drop"
  )

rawCount_wide_combined <- rawCount_stats %>%
  mutate(
    mean_r  = round(mean_count, 2),
    se_r    = ifelse(is.na(se_count), NA, round(se_count, 2)),
    mean_se = ifelse(is.na(se_r), as.character(mean_r), paste0(mean_r, " \u00b1 ", se_r))
  ) %>%
  select(name, Treatment, mean_se) %>%
  pivot_wider(names_from = Treatment, values_from = mean_se, values_fill = "") %>%
  distinct(name, .keep_all = TRUE) %>%
  column_to_rownames("name")

rawCount_wide_numeric <- rawCount_stats %>%
  select(name, Treatment, mean_count, se_count) %>%
  pivot_wider(names_from = Treatment,
              values_from = c(mean_count, se_count),
              names_sep = "_", values_fill = NA) %>%
  distinct(name, .keep_all = TRUE) %>%
  column_to_rownames("name")

rawCount_combined_df <- rawCount_wide_combined %>% rownames_to_column(var = "GENE")
rawCount_numeric_df  <- rawCount_wide_numeric  %>% rownames_to_column(var = "GENE")

merged_all2 <- merged_all %>%
  left_join(rawCount_combined_df, by = "GENE") %>%
  left_join(rawCount_numeric_df,  by = "GENE")

out_final <- paste0(out_prefix, "_merged_pairwise_with_overall_and_rawcounts.tsv")
write.table(merged_all2, file = out_final, sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved merged table with raw-count means and SE to: ", out_final)

added_cols <- setdiff(colnames(merged_all2), colnames(merged_all))
message("Columns added: ", paste(head(added_cols, 20), collapse = ", "))
print(head(merged_all2[, c("GENE", added_cols)], 6))

# ---------- Pairwise plots ----------
out_pdf          <- "output/aldex_pairwise_plots_and_GENE_plots.pdf"
out_csv_prefix   <- "output/aldex_pairwise_publication_table"
effect_threshold <- 0.8
overlap_threshold <- 0.20
top_n_gene_plots  <- 5

stopifnot(exists("counts_f"), exists("meta"))
stopifnot(identical(colnames(counts_f), rownames(meta)))
if (!exists("comparisons") || length(comparisons) == 0) {
  stop("Please ensure 'comparisons' list exists.")
}

raw_results <- list()
for (comp in names(comparisons)) {
  groups_vec  <- comparisons[[comp]]
  keep_samples <- rownames(meta)[meta$Treatment %in% groups_vec]
  meta_sub    <- meta[keep_samples, , drop = FALSE]
  counts_sub  <- counts_f[, keep_samples, drop = FALSE]
  message("Preparing ALDEx2 object for ", comp, " (", length(keep_samples), " samples)")
  raw_results[[comp]] <- aldex(reads      = counts_sub,
                               conditions = meta_sub$Treatment,
                               mc.samples = mc.samples_pairwise,
                               test       = "t",
                               effect     = TRUE,
                               CI         = FALSE,
                               include.sample.summary = FALSE,
                               verbose    = FALSE,
                               paired.test = FALSE,
                               denom      = denom_choice,
                               iterate    = FALSE,
                               gamma      = NULL)
}

pdf(out_pdf, width = 14, height = 8)
message("Writing plots to PDF: ", out_pdf)

for (comp in names(raw_results)) {
  res_obj <- raw_results[[comp]]
  res_df  <- as.data.frame(res_obj)
  res_df$GENE <- rownames(res_df)

  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  tryCatch(
    aldex.plot(res_obj, type = "MA", test = "welch",
               xlab = "Log-ratio abundance", ylab = "Difference",
               main = paste0(comp, " \u2014 Bland-Altman (MA)")),
    error = function(e) { plot.new(); title(main = paste0("MA plot failed: ", comp)) }
  )
  tryCatch(
    aldex.plot(res_obj, type = "MW", test = "welch",
               xlab = "Dispersion", ylab = "Difference",
               main = paste0(comp, " \u2014 Effect (MW)")),
    error = function(e) { plot.new(); title(main = paste0("MW plot failed: ", comp)) }
  )
  tryCatch(
    aldex.plot(res_obj, type = "volcano", test = "welch",
               xlab = "Difference", ylab = "-1(log10(q))",
               main = paste0(comp, " \u2014 Volcano")),
    error = function(e) { plot.new(); title(main = paste0("Volcano plot failed: ", comp)) }
  )

  qcol <- if ("wi.eBH" %in% colnames(res_df)) "wi.eBH" else if ("we.eBH" %in% colnames(res_df)) "we.eBH" else NULL
  if (is.null(qcol)) {
    plot.new(); title(main = paste0(comp, " \u2014 No BH q column found"))
  } else {
    plot(res_df$effect, -log10(res_df[[qcol]]), pch = 20,
         xlab = "ALDEx2 effect (CLR diff)", ylab = paste0("-log10(", qcol, ")"),
         main = paste0(comp, " \u2014 Volcano (effect vs -log10 q)"))
    abline(v = c(-effect_threshold, effect_threshold), col = "grey", lty = 2)
    abline(h = -log10(0.05), col = "red", lty = 2)
  }

  cand   <- res_df %>% filter(!is.na(effect)) %>%
    filter(abs(effect) >= effect_threshold & overlap <= overlap_threshold) %>%
    arrange(desc(abs(effect)))

  qcol_w <- if ("wi.eBH" %in% colnames(cand)) "wi.eBH" else if ("we.eBH" %in% colnames(cand)) "we.eBH" else NA
  pub_cols <- unique(c("GENE", "effect", "overlap", qcol_w,
                        grep("^rab\\.win\\.", colnames(cand), value = TRUE)))
  pub_df <- cand %>%
    select(any_of(pub_cols)) %>%
    mutate(effect = round(effect, 3), overlap = round(overlap, 4),
           across(all_of(qcol_w), ~ round(., 4)))
  out_csv <- paste0(out_csv_prefix, "_", comp, ".csv")
  write_csv(pub_df, out_csv)
  message("Wrote publication CSV for ", comp, ": ", out_csv)

  if (nrow(cand) == 0) {
    top_gene <- res_df %>% arrange(desc(abs(effect))) %>% slice_head(n = top_n_gene_plots) %>% pull(GENE)
    message("No candidates passed thresholds for ", comp, "; plotting top ", length(top_gene), " gene by effect.")
  } else {
    top_gene <- head(cand$GENE, top_n_gene_plots)
    message("Plotting ", length(top_gene), " candidate gene for ", comp)
  }

  groups_vec   <- comparisons[[comp]]
  keep_samples <- rownames(meta)[meta$Treatment %in% groups_vec]
  meta_sub     <- meta[keep_samples, , drop = FALSE]
  counts_sub   <- counts_f[, keep_samples, drop = FALSE]

  counts_long <- as.data.frame(counts_sub) %>%
    rownames_to_column(var = "GENE") %>%
    pivot_longer(cols = -GENE, names_to = "sample", values_to = "count") %>%
    left_join(meta_sub %>% rownames_to_column(var = "sample"), by = "sample") %>%
    mutate(log_count = log10(count + 1))

  for (gen in top_gene) {
    class_data <- counts_long %>% filter(GENE == gen)
    if (nrow(class_data) == 0) next

    kw_res       <- kruskal.test(count ~ Treatment, data = class_data)
    sample_stats <- class_data %>%
      group_by(sample, Treatment) %>%
      summarise(mean_count = mean(count, na.rm = TRUE),
                se_count   = ifelse(n() > 1, sd(count, na.rm = TRUE) / sqrt(n()), NA_real_),
                .groups = "drop")

    barp <- ggplot(class_data, aes(x = sample, y = count, fill = Treatment)) +
      geom_col() +
      geom_text(aes(label = round(count, 1)), vjust = -1.6, size = 3) +
      geom_hline(aes(yintercept = mean(count)), linetype = "dashed", color = "red") +
      annotate("text", x = Inf, y = mean(class_data$count), label = "",
               hjust = 1.1, vjust = -0.5, color = "red", size = 3) +
      geom_line(aes(group = 1), color = "black", linewidth = 1) +
      geom_errorbar(data = sample_stats,
                    aes(x = sample, ymin = pmax(0, mean_count - se_count), ymax = mean_count + se_count),
                    width = 0.2, color = "black", inherit.aes = FALSE) +
      geom_point(data = sample_stats, aes(x = sample, y = mean_count),
                 color = "red", size = 2, inherit.aes = FALSE) +
      scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                                   "Dulce"          = "#2ca02c",
                                   "Soyabean meal"  = "#ff7f0e")) +
      facet_wrap(~ Treatment, scales = "free_x", nrow = 1) +
      labs(title = bquote(italic(.(gen)) ~ "counts per sample (KW p =" ~
                            .(signif(kw_res$p.value, 3)) ~ ")"),
           x = "Sample", y = "Raw count") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title  = element_text(hjust = 0.5))

    present_groups <- unique(class_data$Treatment)
    my_comparisons <- list()
    if (all(c("Reference diet", "Dulce") %in% present_groups))
      my_comparisons <- append(my_comparisons, list(c("Reference diet", "Dulce")))
    if (all(c("Soyabean meal", "Reference diet") %in% present_groups))
      my_comparisons <- append(my_comparisons, list(c("Soyabean meal", "Reference diet")))
    if (all(c("Dulce", "Soyabean meal") %in% present_groups))
      my_comparisons <- append(my_comparisons, list(c("Dulce", "Soyabean meal")))

    label_KW <- paste0("Kruskal-Wallis p = ", signif(kw_res$p.value, 3))

    violp <- ggplot(class_data, aes(x = Treatment, y = log_count, fill = Treatment)) +
      geom_violin(trim = FALSE, scale = "width") +
      geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
      scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                                   "Dulce"          = "#2ca02c",
                                   "Soyabean meal"  = "#ff7f0e")) +
      labs(x = "Treatment", y = "log10(count)") +
      theme_classic() +
      theme(legend.position = "none", text = element_text(size = 12)) +
      annotate("text", x = 1.5, y = max(class_data$log_count, na.rm = TRUE) + 0.5,
               label = label_KW, size = 4)

    if (length(my_comparisons) > 0) {
      y_max   <- max(class_data$log_count, na.rm = TRUE)
      label_y <- seq(y_max + 0.5, y_max + 0.5 + 0.6 * (length(my_comparisons) - 1), by = 0.6)
      violp   <- violp + stat_compare_means(comparisons = my_comparisons,
                                            method = "wilcox.test",
                                            label  = "p.signif",
                                            label.y = label_y)
    }
    print(barp)
    print(violp)
  }
}

dev.off()
message("All plots written to: ", out_pdf)

# -------------------- Top gene summary --------------------
out_pdf2         <- "output/aldex_top_gene_bar_violin.pdf"
out_csv_summary  <- "output/aldex_top_candidates_summary.csv"
effect_threshold  <- 0.8
overlap_threshold <- 0.20
per_comp_candidates <- 20
plot_top_n          <- 10
pairwise_comps <- c("Dulce_vs_Reference", "Reference_vs_Soybean", "Dulce_vs_Soybean")

stopifnot(exists("merged_all2"), exists("counts_f"), exists("meta"))
counts_mat <- counts_f
stopifnot(identical(colnames(counts_mat), rownames(meta)))

get_comp_cols <- function(merged_df, comp_prefix) {
  cols        <- colnames(merged_df)
  effect_col  <- grep(paste0("^", comp_prefix, "_.*effect$"), cols, value = TRUE)
  overlap_col <- grep(paste0("^", comp_prefix, "_.*overlap$"), cols, value = TRUE)
  q_candidates <- c(paste0(comp_prefix, "_wi.eBH"), paste0(comp_prefix, "_we.eBH"),
                    paste0(comp_prefix, "_wi.ep"),  paste0(comp_prefix, "_we.ep"))
  q_col <- intersect(q_candidates, cols)
  if (length(q_col) == 0)
    q_col <- grep(paste0("^", comp_prefix, ".*(eBH|ep)$"), cols, value = TRUE)
  q_col <- if (length(q_col) > 0) q_col[1] else NA_character_
  list(effect  = ifelse(length(effect_col) > 0,  effect_col[1],  NA_character_),
       overlap = ifelse(length(overlap_col) > 0, overlap_col[1], NA_character_),
       q       = q_col)
}

top_candidates_list <- list()
for (comp in pairwise_comps) {
  cols <- get_comp_cols(merged_all2, comp)
  if (is.na(cols$effect) || is.na(cols$overlap)) {
    message("Warning: could not find effect/overlap columns for ", comp, ". Skipping.")
    next
  }
  tmp <- merged_all2 %>%
    select(GENE, all_of(c(cols$effect, cols$overlap, cols$q))) %>%
    rename(effect = all_of(cols$effect), overlap = all_of(cols$overlap))
  if (!is.na(cols$q)) names(tmp)[names(tmp) == cols$q] <- "pair_q"
  tmp  <- tmp %>% mutate(across(c(effect, overlap, pair_q), as.numeric))
  cand <- tmp %>% filter(!is.na(effect)) %>%
    filter(abs(effect) >= effect_threshold & overlap <= overlap_threshold) %>%
    arrange(desc(abs(effect)))
  if (nrow(cand) < per_comp_candidates) {
    more <- tmp %>% filter(!is.na(effect)) %>%
      arrange(desc(abs(effect))) %>% slice_head(n = per_comp_candidates)
    cand <- distinct(bind_rows(cand, more), GENE, .keep_all = TRUE)
  } else {
    cand <- cand %>% slice_head(n = per_comp_candidates)
  }
  top_candidates_list[[comp]] <- cand$GENE
  message("Selected ", length(cand$GENE), " candidates for ", comp)
}

all_candidates <- unique(unlist(top_candidates_list))
message("Total unique candidates across comparisons: ", length(all_candidates))

if (length(all_candidates) > (plot_top_n * length(pairwise_comps))) {
  effect_cols_all <- unlist(lapply(pairwise_comps, function(p) get_comp_cols(merged_all2, p)$effect))
  effect_cols_all <- effect_cols_all[!is.na(effect_cols_all)]
  rank_df <- merged_all2 %>%
    filter(GENE %in% all_candidates) %>%
    mutate(max_abs_effect = apply(abs(select(., any_of(effect_cols_all))), 1, max, na.rm = TRUE)) %>%
    arrange(desc(max_abs_effect))
  keep_n         <- min(length(all_candidates), plot_top_n * length(pairwise_comps))
  all_candidates <- rank_df$GENE[1:keep_n]
  message("Reduced candidates to top ", keep_n, " by max abs(effect).")
}

counts_long <- as.data.frame(counts_mat) %>%
  rownames_to_column(var = "GENE") %>%
  pivot_longer(cols = -GENE, names_to = "sample", values_to = "count") %>%
  left_join(meta %>% rownames_to_column(var = "sample"), by = "sample") %>%
  mutate(log_count = log10(count + 1))

pdf(out_pdf2, width = 12, height = 8)
message("Writing plots to: ", out_pdf2)

for (gen in all_candidates) {
  class_data <- counts_long %>% filter(GENE == gen)
  if (nrow(class_data) == 0) next

  kw_res   <- tryCatch(kruskal.test(count ~ Treatment, data = class_data), error = function(e) NULL)
  kw_label <- if (!is.null(kw_res)) paste0("Kruskal-Wallis p = ", signif(kw_res$p.value, 3)) else "KW test failed"

  sample_stats <- class_data %>%
    group_by(sample, Treatment) %>%
    summarise(mean_count = mean(count, na.rm = TRUE),
              se_count   = ifelse(n() > 1, sd(count, na.rm = TRUE) / sqrt(n()), NA_real_),
              .groups = "drop")

  barp <- ggplot(class_data, aes(x = reorder(sample, count), y = count, fill = Treatment)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(count, 1)), vjust = -1.6, size = 3) +
    geom_hline(aes(yintercept = mean(count)), linetype = "dashed", color = "red") +
    annotate("text", x = Inf, y = mean(class_data$count), label = "",
             hjust = 1.1, vjust = -0.5, color = "red", size = 3) +
    geom_line(aes(group = 1), color = "black", linewidth = 1) +
    geom_errorbar(data = sample_stats,
                  aes(x = sample, ymin = mean_count - se_count, ymax = mean_count + se_count),
                  width = 0.2, color = "white", inherit.aes = FALSE) +
    geom_point(data = sample_stats, aes(x = sample, y = mean_count),
               color = "red", size = 2, inherit.aes = FALSE) +
    scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                                 "Dulce"          = "#2ca02c",
                                 "Soyabean meal"  = "#ff7f0e")) +
    labs(title = bquote(italic(.(gen)) ~ "Count per Treatment (" ~ .(kw_label) ~ ")"),
         x = "Sample groups", y = "Count") +
    theme_classic() +
    theme(legend.position = "none", text = element_text(size = 14),
          axis.text.y    = element_text(angle = 90, hjust = 0.5),
          plot.title     = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x    = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~ Treatment, scales = "free_x", nrow = 1)

  my_comparisons <- list(c("Reference diet", "Dulce"),
                         c("Soyabean meal",  "Reference diet"),
                         c("Dulce",          "Soyabean meal"))
  present_groups <- unique(class_data$Treatment)
  my_comparisons <- Filter(function(x) all(x %in% present_groups), my_comparisons)

  violp <- ggplot(class_data, aes(x = Treatment, y = log_count, fill = Treatment)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
    scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                                 "Dulce"          = "#2ca02c",
                                 "Soyabean meal"  = "#ff7f0e")) +
    labs(x = "Treatment", y = "log10(count + 1)") +
    theme_classic() +
    theme(legend.position = "none", text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0(gen, " \u2014 group comparison")) +
    annotate("text", x = 2, y = max(class_data$log_count, na.rm = TRUE) + 1.0,
             label = kw_label, size = 4)

  if (length(my_comparisons) > 0) {
    y_max   <- max(class_data$log_count, na.rm = TRUE)
    label_y <- seq(y_max + 0.2, y_max + 0.2 + 0.3 * (length(my_comparisons) - 1), by = 0.3)
    violp   <- violp + stat_compare_means(comparisons = my_comparisons,
                                          method  = "wilcox.test",
                                          label   = "p.signif",
                                          label.y = label_y)
  }
  print(barp)
  print(violp)
}

dev.off()
message("PDF saved: ", out_pdf2)

# Publication-ready CSV for selected genes
overall_cols <- grep("kw\\.|kw\\b|kw\\.eBH|kw\\.ep", colnames(merged_all2), value = TRUE)
if (length(overall_cols) == 0)
  overall_cols <- grep("kw|KW|overall", colnames(merged_all2), value = TRUE)

pairwise_summary_cols <- unlist(lapply(pairwise_comps, function(p) {
  cols <- get_comp_cols(merged_all2, p)
  c(cols$effect, cols$overlap, cols$q)
}))
pairwise_summary_cols <- pairwise_summary_cols[!is.na(pairwise_summary_cols)]
summary_cols <- unique(c("GENE", overall_cols, pairwise_summary_cols))

pub_table <- merged_all2 %>%
  filter(GENE %in% all_candidates) %>%
  select(any_of(summary_cols)) %>%
  mutate(across(where(is.numeric), ~ round(., 4)))

write_csv(pub_table, out_csv_summary)
message("Publication CSV written: ", out_csv_summary)

message("All done.")

# chickMicro R Scripts — Pre-GitHub Cleanup Log

**Date:** 2026-03-30
**Scope:** `chickMicro/R_scripts/`
**Reason:** Remove non-project content, fix hardcoded paths, remove dead code, and fix bugs before public GitHub submission.

---

## Summary Table

| File | Action | Details |
|------|--------|---------|
| `one-way_ANOVA.R` | DELETED | Used simulated/example data unrelated to the project |
| `kraken2_VennDiagram.R` | REWRITTEN | Was a Lung microbiome script (Exhale vs Sputum); replaced with chicken microbiome version comparing 3 treatment groups using `abri_kraken2_filtered.csv` |
| `relativeAbundance.R` | REWRITTEN | Was a Lung microbiome script referencing `noncontaminants_list.csv`; replaced with chicken microbiome stacked bar chart using `abri_kraken2_filtered.csv` and `chicken_metadata.csv` |
| `cleanData.R` | PATH FIXES + DEAD CODE REMOVED | 6 absolute paths → relative; large commented block (lines 91–190) removed |
| `normData.R` | PATH FIXES + DEAD CODE REMOVED | 2 absolute paths → relative; large commented pipeline block (lines 8–41) removed |
| `chordDiagram_v1.2.1.R` | PATH FIXES | 2 absolute paths → relative |
| `clostridiaCount_Treatment.R` | PATH FIXES | 2 absolute paths → relative |
| `geneMap.R` | PATH FIXES | 1 absolute path → relative |
| `richnessTreatment.R` | PATH FIXES | 1 absolute path → relative |
| `statsAnova.R` | PATH FIXES + BUG FIX | 2 absolute paths → relative; merge() moved before the filter block |
| `taxaNormProfile_PCoA.R` | PATH FIXES | 2 absolute paths → relative |
| `taxaNorm_pheatmap.R` | PATH FIXES + DEAD CODE REMOVED | 3 absolute paths → relative; commented `rearrange_matrix_by_category` block removed |
| `taxaProfile_PCoA.R` | PATH FIXES + COLUMN INDEX FIX | 4 absolute paths → relative (including saveWidget and write.csv); hardcoded `2:19` column indices replaced with `everything()` / `-name` |
| `tripartiteNetwork.R` | PATH FIXES | 3 absolute paths → relative |
| `violinNormPlot.R` | PATH FIXES | 1 absolute path → relative |
| `violinPlot.R` | PATH FIXES + DEAD CODE REMOVED + BUG FIX | 2 absolute paths → relative; commented bracken section removed; dead AMR block using undefined `abri_kraken2Bracken_merged` removed |

---

## Section 1: Deleted Files

### `one-way_ANOVA.R`
- **Reason:** This script used simulated/example data and was not connected to any project datasets. It had no scientific contribution to the chicken gut microbiome analysis.
- **Action:** Permanently deleted.

---

## Section 2: Rewritten Files

### `kraken2_VennDiagram.R`
- **Problem:** The file was accidentally copied from a Lung microbiome project. It loaded `noncontaminants_list.csv` from the Lung microbiome Re-centrifuge directory, compared "Exhale" vs "Sputum" sample groups, and used lowercase `taxid` as the identifier.
- **Replacement:** New file loads `../data/abri_kraken2_filtered.csv`, splits data into three treatment groups ("Reference diet", "Seaweed", "Soyabean meal"), uses uppercase `TAXID` as the taxon identifier and `TREATMENT` for group splitting. Modelled after `abri_kraken2_VennDiagram.R`. A header comment was added explaining purpose, input, and output.

### `relativeAbundance.R`
- **Problem:** The file was accidentally copied from a Lung microbiome project. It loaded `noncontaminants_list.csv` and `Bracken_noncontaminants_list.csv` from the Lung microbiome directory, applied manual row edits for unrelated species, used "Exhale Breath and Sputum Sample Groups" as the x-axis label, and saved output to a Lung microbiome figures directory.
- **Replacement:** New file loads `../data/abri_kraken2_filtered.csv` and `../data/chicken_metadata.csv`, merges on `SAMPLE`, filters and abbreviates species names consistent with other project scripts, computes per-sample relative abundance percentages, collapses to top-50 taxa, uses the `paletteer` palette approach, and renders a stacked bar chart faceted by `TREATMENT`. A header comment was added explaining purpose, input, and output.

---

## Section 3: Path Fixes (Absolute → Relative)

All scripts run from `R_scripts/`, so data files are addressed via `"../data/"` and figure outputs via `"../figures/"`.

### `cleanData.R`
| Old path | New path |
|----------|----------|
| `C:/.../Reports/Abricate/summary_reportClean.csv` | `"../data/summary_reportClean.csv"` |
| `C:/.../Reports/Kraken2/kraken2_db1_combined_reports.txt` | `"../data/kraken2_db1_combined_reports.txt"` |
| `C:/.../Reports/Kraken2/kraken2_db2-bracken_combined_reports.txt` | `"../data/kraken2_db2-bracken_combined_reports.txt"` |
| `C:/.../Reports/Re-centrifuge/PC_JC_2024-11-29_Julio_Gallus.rcf.data.csv` | `"../data/recentrifuge_contaminants.csv"` |
| `C:/.../Datasets/abri_kraken2_cleaned.csv` (write output) | `"../data/abri_kraken2_filtered.csv"` |
| `C:/.../Reports/Kraken2/bracken_arranged.csv` (write output) | `"../data/bracken_arranged.csv"` |

### `normData.R`
| Old path | New path |
|----------|----------|
| `C:/.../Datasets/abri_kraken2Bracken_merged.csv` | `"../data/abri_kraken2Bracken_merged.csv"` |
| `C:/.../Datasets/genetable_normdata.csv` (write output) | `"../data/genetable_normdata.csv"` |

### `chordDiagram_v1.2.1.R`
| Old path | New path |
|----------|----------|
| `C:/.../Datasets/abri_kraken2_cleaned.csv` | `"../data/abri_kraken2_filtered.csv"` |
| `C:/.../Datasets/chicken_metadata.csv` | `"../data/chicken_metadata.csv"` |

### `clostridiaCount_Treatment.R`
| Old path | New path |
|----------|----------|
| `C:/.../Reports/Kraken2/bracken_arranged.csv` | `"../data/bracken_arranged.csv"` |
| `C:/.../Datasets/chicken_metadata.csv` | `"../data/chicken_metadata.csv"` |

### `geneMap.R`
| Old path | New path |
|----------|----------|
| `C:/.../Datasets/abri_kraken2_merged.csv` | `"../data/abri_kraken2_merged.csv"` |

### `richnessTreatment.R`
| Old path | New path |
|----------|----------|
| `C:/.../Datasets/genetable_normdata.csv` | `"../data/genetable_normdata.csv"` |

### `statsAnova.R`
| Old path | New path |
|----------|----------|
| `C:/.../Datasets/abri_kraken2_cleaned.csv` | `"../data/abri_kraken2_filtered.csv"` |
| `C:/.../Datasets/chicken_metadata.csv` | `"../data/chicken_metadata.csv"` |

### `taxaNormProfile_PCoA.R`
| Old path | New path |
|----------|----------|
| `C:/.../Reports/Kraken2/bracken_arranged.csv` | `"../data/bracken_arranged.csv"` |
| `C:/.../Datasets/chicken_metadata.csv` | `"../data/chicken_metadata.csv"` |

### `taxaNorm_pheatmap.R`
| Old path | New path |
|----------|----------|
| `C:/.../Reports/Kraken2/bracken_arranged.csv` | `"../data/bracken_arranged.csv"` |
| `C:/.../Datasets/chicken_metadata.csv` | `"../data/chicken_metadata.csv"` |
| `C:/.../aldex_top_candidates_summary.csv` | `"../data/aldex_top_candidates_summary.csv"` |

### `taxaProfile_PCoA.R`
| Old path | New path |
|----------|----------|
| `C:/.../Reports/Kraken2/bracken_arranged.csv` | `"../data/bracken_arranged.csv"` |
| `C:/.../Datasets/chicken_metadata.csv` | `"../data/chicken_metadata.csv"` |
| `C:/.../Figures/PCoA_plot.html` (saveWidget output) | `"../figures/PCoA_plot.html"` |
| `C:/.../Datasets/ann_pcoa.csv` (write.csv output) | `"../data/ann_pcoa.csv"` |

### `tripartiteNetwork.R`
| Old path | New path |
|----------|----------|
| `C:/.../Datasets/abri_kraken2Bracken_merged.csv` | `"../data/abri_kraken2Bracken_merged.csv"` |
| `C:/.../Datasets/gephi_edges3.csv` (write output) | `"../data/gephi_edges3.csv"` |
| `C:/.../Datasets/gephi_nodes3.csv` (write output) | `"../data/gephi_nodes3.csv"` |

### `violinNormPlot.R`
| Old path | New path |
|----------|----------|
| `C:/.../Datasets/genetable_normdata.csv` | `"../data/genetable_normdata.csv"` |

### `violinPlot.R`
| Old path | New path |
|----------|----------|
| `C:/.../Datasets/abri_kraken2_cleaned.csv` | `"../data/abri_kraken2_filtered.csv"` |
| `C:/.../Datasets/chicken_metadata.csv` | `"../data/chicken_metadata.csv"` |

---

## Section 4: Dead Commented Code Removed

### `cleanData.R`
- **Removed:** Large block from `# # Clean up contaminants data set` through `# write.csv(abri_kraken2_uniqued, ...)` (original lines ~91–190). This block contained an abandoned workflow for processing Re-centrifuge contaminant/non-contaminant data that was superseded by the `abri_kraken2_filtered.csv` pipeline. It included multiple `pivot_longer`, `filter`, and `write.csv` operations referencing no-longer-used intermediate files.

### `normData.R`
- **Removed:** Entire commented block at top of file (original lines ~8–46). This block contained an earlier version of the data loading and merge pipeline (reading `abri_kraken2_merged.csv`, `bracken_arranged.csv`, performing pivot and merge to create `abri_kraken2Bracken_merged`). It was superseded by loading the pre-built `abri_kraken2Bracken_merged.csv` directly (which remains on line 49).

### `taxaNorm_pheatmap.R`
- **Removed:** Commented-out `rearrange_matrix_by_category` function definition and its call (`reordered_matrix <- rearrange_matrix_by_category(...)`). This function was replaced by the `rearrange_matrix` function used elsewhere in the same script, and the reordered matrix variable was never used downstream.

### `violinPlot.R`
- **Removed (commented bracken import):** Three commented lines loading `bracken_arranged` from absolute path (`# bracken_arranged <- read_csv(...)`, `# View(...)`, `# colnames(...)`).
- **Removed (commented bracken pivot block):** Seven commented lines for `bracken_pivoted` and `bracken_merged` that were part of an abandoned workflow using the Bracken report instead of the Abricate output.
- **Removed (dead AMR block):** The entire block from `# Filter the data set by AMR` through `coord_flip()` (original lines ~100–132). This block referenced the undefined variable `abri_kraken2Bracken_merged`, which was never loaded in this script. The block computed `Gene_log_count` from a non-existent `Gene_count` column and would error at runtime.

---

## Section 5: Bug Fixes

### `statsAnova.R` — Undefined variable (`abri_kraken2_merged` used before creation)
- **Problem:** Line ~39 assigned `abri_kraken2_filtered <- abri_kraken2_merged %>% filter(...)` but `abri_kraken2_merged` was not created until line ~57 (`abri_kraken2_merged <- merge(abri_kraken2_filtered, chicken_metadata, by = "sample")`). Running the script would fail with "object 'abri_kraken2_merged' not found".
- **Fix:** Moved the `merge()` call to immediately after the metadata string-cleaning step (before the filter block). The merge now uses the originally loaded `abri_kraken2_filtered` (from `read_csv`) and `chicken_metadata`. The filter block then correctly operates on `abri_kraken2_merged`.

### `taxaProfile_PCoA.R` — Hardcoded column indices
- **Problem:** Lines ~68–69 used `across(2:19, ...)` in both `summarise()` and `mutate()` calls. These indices assumed exactly 18 sample columns (D19–D36), making the script fragile and liable to fail or produce wrong results if the number of samples changes.
- **Fix:**
  - `summarise(across(2:19, ...))` → `summarise(across(everything(), ...))`
  - `mutate(across(2:19, ...))` → `mutate(across(-name, ...))`
  - The grouping variable in this script is `name` (lowercase), so `-name` correctly excludes the identifier column while applying the string-replace operation to all sample columns dynamically.

### `violinPlot.R` — Undefined variable (`abri_kraken2Bracken_merged`)
- **Problem:** Line ~101 used `abri_kraken2Bracken_merged` which was never loaded in `violinPlot.R`. This would cause an immediate "object not found" error.
- **Fix:** The entire dead AMR block that referenced this variable was removed (see Section 4). The rest of the script that uses the correctly-loaded `abri_kraken2_merged` (from the merge of `abri_kraken2_clean` and `chicken_metadata`) was preserved intact.

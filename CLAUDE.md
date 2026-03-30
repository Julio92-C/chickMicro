# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Is

**chickMicro** is an R Shiny web application for analyzing poultry gut microbiome metagenomic data, with a focus on antimicrobial resistance (AMR) profiling. It processes outputs from Abricate (gene detection) and Kraken2/Bracken (taxonomic classification).

Live deployment: https://julio92-c.shinyapps.io/hosMicro/

## Running the App

```r
# From R or RStudio, set working directory to the app folder:
setwd("chickMicro/app")
shiny::runApp()
```

Or open `chickMicro/app/ui.R` in RStudio and click "Run App".

## Deploying to shinyapps.io

```r
library(rsconnect)
rsconnect::deployApp("chickMicro/app", appName = "hosMicro")
```

## Running Standalone Analysis Scripts

```r
# From within R:
source("chickMicro/R_scripts/cleanData.R")

# Or from terminal:
Rscript chickMicro/R_scripts/cleanData.R
```

## Architecture

The app uses the **Shiny module pattern** for separation of concerns. Entry points are:

- `app/global.R` — sources all modules on startup
- `app/ui.R` — dashboard layout with 4 tabs (Home, Dashboard, Contact, About)
- `app/server.R` — authentication layer + wires module servers together

### Module System (`app/modules/`)

| Module file | UI function | Server function | Purpose |
|---|---|---|---|
| `input_csvFile_Module.R` | `csvFileUI()` | `csvFileServer()` | CSV upload and parsing |
| `mdr_dataPlot_Module.R` | `MDRdataPlotUI()` | `MDRdataPlotServer()` | Filtering, visualization |
| `dataStats_Module.R` | `dataStatsUI()` | `dataStatsServer()` | Summary statistics tables |

### Reactive Data Flow

```
User uploads CSV → csvFileServer (reactive dataframe)
    → MDRdataPlotServer (database filter + threshold filters → plots)
    → dataStatsServer (filtered data → gt summary table)
```

Filters available in MDRdataPlotServer: database (CARD/VFDB/PlasmidFinder), min taxa count, min gene count, coverage %, identity %.

### Report Generation

`app/reports/shinyReport.Rmd` — R Markdown with `runtime: shiny`, receives the filtered dataframe as a parameter and generates an HTML report with interactive DataTables.

## Key Data Structures

Input CSVs are expected to have these columns (from Abricate + Kraken2 pipeline):
`TAXID, SAMPLE, SEQUENCE_ID, GENE, COVERAGE_PCT, IDENTITY_PCT, DATABASE, RESISTANCE, NAME, PRODUCT`

Metadata CSV columns: `SAMPLE, TREATMENT`

Sample treatments in the study: `"Reference diet"`, `"Soyabean meal"`, `"Seaweed"`

## Authentication

`server.R` has hardcoded demo credentials (`admin/admin123`, `guest/guest123`). This is intentional for the demo deployment — do not replace with a production auth system unless the deployment context changes.

## R Scripts (Non-Shiny Analysis)

`R_scripts/` contains 30+ standalone scripts for the underlying research pipeline:
- **Data prep:** `cleanData.R`, `normData.R`
- **Diversity:** `ARGNorm_alphaDiversity.R`, `ARGdiversity_check.R`
- **Visualization:** VennDiagram scripts, pheatmap scripts, PCoA, Sankey diagrams
- **Stats:** `statsAnova.R`, `one-way_ANOVA.R`, `AMR_VFs_MGEs_Stats_Summary.R`

These are research scripts, not part of the Shiny app runtime.

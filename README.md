# 🐔 chickMicro

> An R Shiny web application to explore the chicken gut microbiome and their genetic makeup.

---

## 🔬 About

**chickMicro** is an interactive web application for analysing poultry gut microbiome metagenomic data, with a focus on antimicrobial resistance (AMR) profiling. It processes outputs from **Abricate** (gene detection) and **Kraken2/Bracken** (taxonomic classification).

### ✨ Features

- 📂 Upload and parse CSV files generated from metagenomic pipelines
- 🔍 Filter data by database (CARD, VFDB, PlasmidFinder), coverage %, identity %, and gene/taxa count thresholds
- 📊 Visualise AMR gene profiles and multidrug resistance (MDR) patterns across samples
- 📋 Explore summary statistics in interactive tables
- 📄 Generate downloadable HTML reports from filtered data

### 🧪 Dataset

The dataset comprises broiler chicken gut samples across three dietary treatment groups: **Reference diet**, **Soyabean meal**, and **Seaweed**.

---

## 🚀 Running the App

```r
# From R or RStudio, set working directory to the app folder:
setwd("chickMicro/app")
shiny::runApp()
```

Or open `chickMicro/app/ui.R` in RStudio and click **Run App**.

---

## 📬 Contact

👤 **Julio C. Ortega Cambara**

🎓 PhD Candidate — Computational Bioinformatics

🏛️ School of Biomedical Sciences, University of West London

📍 St Mary's Rd, London W5 5RF

✉️ Email: [32104617@student.uwl.ac.uk](mailto:32104617@student.uwl.ac.uk)

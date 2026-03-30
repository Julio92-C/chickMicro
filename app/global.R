
source("modules/input_csvFile_Module.R")
source("modules/mdr_dataPlot_Module.R")
source("modules/dataStats_Module.R")

# Pre-load datasets on app startup
app_data     <- read.csv("data/abri_kraken2_filtered.csv", stringsAsFactors = FALSE)
app_metadata <- read.csv("data/chicken_metadata.csv",      stringsAsFactors = FALSE)

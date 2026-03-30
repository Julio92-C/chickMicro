# install_packages.R
# Run this script once to install all packages required by chickMicro.
# In RStudio: open this file and click "Source", or run:
#   source("chickMicro/install_packages.R")

pkgs <- c(
  # Shiny & dashboard
  "shiny",
  "shinydashboard",
  "shinydashboardPlus",
  "shinythemes",
  "shinycssloaders",
  "thematic",

  # Tables & output
  "DT",
  "gt",
  "gtsummary",

  # Visualisation
  "ggplot2",
  "plotly",
  "pheatmap",
  "VennDiagram",
  "circlize",
  "ggsignif",
  "paletteer",
  "scales",

  # Data wrangling
  "tidyverse",   # includes dplyr, tidyr, stringr, forcats, readr, tibble, purrr
  "dplyr",
  "tidyr",
  "forcats",
  "stringr",
  "reshape2",

  # Statistics / ecology
  "vegan",

  # Reporting
  "rmarkdown"
)

# Install only packages that are not already installed
missing <- pkgs[!pkgs %in% installed.packages()[, "Package"]]

if (length(missing) == 0) {
  message("All packages are already installed.")
} else {
  message("Installing ", length(missing), " missing package(s): ",
          paste(missing, collapse = ", "))
  install.packages(missing, dependencies = TRUE)
  message("Done. Please restart R before running the app.")
}

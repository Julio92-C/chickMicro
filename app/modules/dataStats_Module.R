# Data Stats Module

# Module UI function
dataStatsUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(12,
             box(width = 12, 
                 solidHeader = TRUE, 
                 status = "teal",
                 title = "Basic stats results",
                 collapsible=TRUE,
                 withSpinner(gt_output(ns("my_gt_table")), type = 6, color = "#1a73b8")
             )
             
      )
    )
  )
}




# Module Server function
dataStatsServer <- function(id, filter_data){
  stopifnot(is.reactive(filter_data))
  
  moduleServer(id, function(input, output, session){
    output$my_gt_table <- render_gt({
      validate(
        need(!is.null(filter_data()), "Run a plot first — click the Plot button in the Plot tab."),
        need(nrow(filter_data()) > 0,  "No data to summarise after current filters.")
      )
      # Generate the gt table
      summary_table <- filter_data() %>%
        group_by(SAMPLE) %>%
        summarise(
          DATABASE = list(unique(DATABASE)),
          GENE = list(unique(GENE)),
          COVERAGE_PCT_mean = mean(COVERAGE_PCT),
          COVERAGE_PCT_sd = sd(COVERAGE_PCT),
          IDENTITY_PCT_mean = mean(IDENTITY_PCT),
          IDENTITY_PCT_sd = sd(IDENTITY_PCT)
        )
      summary_table %>%
        gt() %>%
        tab_header(
          title = "Table 1. Summary of Samples groups and genetics elements profile"
        ) %>%
        cols_label(
          DATABASE = "Database",
          GENE = "Gene",
          COVERAGE_PCT_mean = "Mean Coverage (%)",
          COVERAGE_PCT_sd = "Coverage SD (%)",
          IDENTITY_PCT_mean = "Mean Identity (%)",
          IDENTITY_PCT_sd = "Identity SD (%)"
        )
    })
    
  })
  
}
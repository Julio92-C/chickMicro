#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(DT)
library(plotly)
library(ggplot2)
library(dplyr)
library(forcats)
library(thematic)
library(gtsummary)
library(paletteer)
library(scales)
library(reshape2)
library(tidyr)
library(VennDiagram)
library(tidyverse)
library(stringr)
library(ggsignif) 
library(vegan)
library(circlize)
library(pheatmap)

# Login credentials — add or remove users here
credentials <- data.frame(
  user     = c("admin", "guest"),
  password = c("admin123", "guest123"),
  stringsAsFactors = FALSE
)

# Define server logic required to draw a table
shinyServer(function(input, output, session) {
  thematic::thematic_shiny()

  # --- Login state ---
  login_state <- reactiveValues(logged_in = FALSE, username = NULL)

  # Sidebar login/user widget — rendered as native <li> items to match sidebar CSS
  output$loginSidebarUI <- renderUI({
    if (login_state$logged_in) {
      tagList(
        tags$li(
          tags$a(href = "#", icon("user-check"), tags$span(login_state$username))
        ),
        tags$li(
          actionLink("logoutBtn", label = tagList(icon("sign-out-alt"), "Logout"))
        )
      )
    } else {
      tags$li(
        actionLink("loginBtn", label = tagList(icon("user-circle"), "Login"))
      )
    }
  })

  # Show login modal
  observeEvent(input$loginBtn, {
    showModal(modalDialog(
      title = tags$span(icon("user-circle"), " Login"),
      textInput("username_input", "Username", placeholder = "Enter username"),
      passwordInput("password_input", "Password", placeholder = "Enter password"),
      uiOutput("loginError"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("submitLogin", "Login", icon = icon("sign-in-alt"), class = "btn-primary")
      ),
      easyClose = TRUE
    ))
  })

  # Validate credentials
  observeEvent(input$submitLogin, {
    match <- credentials$user == input$username_input &
             credentials$password == input$password_input
    if (any(match)) {
      login_state$logged_in <- TRUE
      login_state$username  <- input$username_input
      removeModal()
    } else {
      output$loginError <- renderUI(
        tags$p(icon("exclamation-triangle"), " Invalid username or password.",
               style = "color: #c0392b; margin-top: 8px;")
      )
    }
  })

  # Logout
  observeEvent(input$logoutBtn, {
    login_state$logged_in <- FALSE
    login_state$username  <- NULL
  })

  # Input csvFileServer Module for MDR-hvKp dataSet
    dataframe <- csvFileServer("MDR_hvKp_dataSet")

    # Input csvFileServer Module for Metadata
    metadata <- csvFileServer("MetadataSet")
  # 
  # # Input MDRdataPlotServer Module for MDR-hvKp dataSet
  filter_data <- MDRdataPlotServer("MDR_hvKp_dataPlot", dataframe = dataframe, metadata = metadata)
  # 
  # # input MDRdataStatsServer Module for MDR-hvKp dataSet
  dataStatsServer("hosMicro_dataStat", filter_data = filter_data)
  
  
  
})

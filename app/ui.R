#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinythemes)
library(plotly)
library(gt)
library(shinydashboard)
library(shinydashboardPlus)
library(rmarkdown)
library(shinycssloaders)

theme = shinythemes::shinytheme("united")  # <--- Specify theme here

# Global variables
source("modules/input_csvFile_Module.R", local = T)
source("modules/mdr_dataPlot_Module.R", local = T)
source("modules/dataStats_Module.R", local = T)



# Define UI for application that draws a histogram
dashboardPage(
  skin = "blue-light",
  dashboardHeader(title = "chickMicro"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home",      tabName = "home",      icon = icon("home")),
      menuItem("Dashboard", tabName = "dashboard",  icon = icon("dashboard")),
      menuItem("Contact",   tabName = "contact",    icon = icon("envelope")),
      menuItem("About",     tabName = "about",      icon = icon("info-circle")),
      uiOutput("loginSidebarUI")
    )
  ),
  
  dashboardBody(
    tags$style("
      /* scrollable body */
      .content-wrapper { overflow-y: auto !important; }

      /* tab pill styling */
      .nav-tabs-custom > .nav-tabs > li > a {
        font-size: 13px;
        font-weight: 600;
        letter-spacing: 0.3px;
        padding: 9px 16px;
        color: #555;
      }
      .nav-tabs-custom > .nav-tabs > li.active > a {
        color: #1a73b8;
        border-top: 3px solid #1a73b8;
      }
      .nav-tabs-custom > .nav-tabs > li > a:hover { color: #1a73b8; }

      /* dashboard header gradient bar */
      .dash-header {
        background: linear-gradient(135deg, #1a73b8, #0d4f8c);
        color: white;
        padding: 16px 22px;
        border-radius: 6px;
        margin-bottom: 14px;
      }
      .dash-header h3 { margin: 0 0 3px 0; font-size: 20px; }
      .dash-header p  { margin: 0; font-size: 12.5px; opacity: 0.85; }

      /* tighten box headers */
      .box.box-solid > .box-header { padding: 8px 14px; }
    "),
    tabItems(
      tabItem(tabName = "home",

        # --- Hero banner ---
        fluidRow(
          column(12,
            tags$div(
              style = "background: linear-gradient(135deg, #1a73b8, #0d4f8c);
                       color: white; padding: 30px 35px; border-radius: 8px; margin-bottom: 20px;",
              tags$h1("Welcome to chickMicro", style = "margin-top: 0;"),
              tags$h4("Exploring the Poultry Gut Microbiome & Antimicrobial Resistance",
                      style = "font-weight: 300; margin-bottom: 12px;"),
              tags$p("A web-based platform for metagenomic data analysis, taxonomy identification,
                      and AMR profiling of poultry production samples.",
                     style = "font-size: 15px; margin-bottom: 0;")
            )
          )
        ),

        # --- Key AMR statistics ---
        fluidRow(
          valueBox(
            value = "10M",
            subtitle = "Projected annual deaths from AMR by 2050",
            icon  = icon("skull-crossbones"),
            color = "red",
            width = 4
          ),
          valueBox(
            value = "1.3M",
            subtitle = "Deaths caused by AMR in 2019",
            icon  = icon("chart-line"),
            color = "orange",
            width = 4
          ),
          valueBox(
            value = "140K",
            subtitle = "Newborn deaths linked to AMR in 2019",
            icon  = icon("baby"),
            color = "yellow",
            width = 4
          )
        ),

        # --- Feature cards + getting started ---
        fluidRow(
          box(
            title = "About This Study", width = 4, status = "primary", solidHeader = TRUE,
            icon  = icon("microscope"),
            tags$p("This application investigates the impact of ",
                   tags$strong("Dulse"), " (a protein-rich red seaweed) as a potential
                   antibiotic alternative on the poultry gut microbiome."),
            tags$p("Understanding these microbial communities is critical for reducing
                   AMR spread in agriculture and protecting public health."),
            tags$p(icon("map-marker-alt"), " University of West London — School of Biomedical Sciences.")
          ),

          box(
            title = "App Features", width = 4, status = "success", solidHeader = TRUE,
            icon  = icon("list"),
            tags$ul(
              style = "padding-left: 18px; line-height: 2;",
              tags$li(icon("upload"),    tags$strong(" Data upload:"),   " CSV metagenomic datasets"),
              tags$li(icon("chart-bar"), tags$strong(" Visualisation:"), " Interactive AMR/VF plots"),
              tags$li(icon("table"),     tags$strong(" Statistics:"),    " Diversity & group analysis"),
              tags$li(icon("circle"),    tags$strong(" Taxonomy:"),      " Gut microbiome profiling"),
              tags$li(icon("dna"),       tags$strong(" MGE profiles:"),  " Mobile genetic elements"),
              tags$li(icon("file-alt"),  tags$strong(" Reports:"),       " Exportable summaries")
            )
          ),

          box(
            title = "Getting Started", width = 4, status = "info", solidHeader = TRUE,
            icon  = icon("play-circle"),
            tags$ol(
              style = "padding-left: 18px; line-height: 2.2;",
              tags$li("Click the ", tags$strong("Dashboard"), " tab in the sidebar"),
              tags$li("Upload your ", tags$strong("CSV data"), " file in the Data panel"),
              tags$li("Upload a ", tags$strong("metadata"), " file if available"),
              tags$li("Switch to ", tags$strong("Plot"), " to explore visualisations"),
              tags$li("Run ", tags$strong("Stats Analysis"), " for diversity metrics")
            )
          )
        ),

        # --- WHO context box ---
        fluidRow(
          box(
            title = "WHO Priority Pathogens & AMR Context", width = 12,
            status = "warning", solidHeader = TRUE, collapsible = TRUE,
            fluidRow(
              column(6,
                tags$h5("Why does this matter?"),
                tags$p("The WHO 2017 Priority Pathogen List identifies bacteria of highest urgency for
                        new antibiotic development. Resistance is driven by overuse in agriculture,
                        inadequate sanitation, and global travel."),
                tags$p("Poultry production represents one of the largest users of antimicrobials
                        worldwide, making the gut microbiome of broiler chickens a key surveillance target.")
              ),
              column(6,
                tags$h5("Research Goals"),
                tags$ul(
                  style = "line-height: 1.9;",
                  tags$li("Characterise the gut microbiome of broiler chickens fed Dulse supplementation"),
                  tags$li("Identify AMR genes, virulence factors (VFs), and mobile genetic elements (MGEs)"),
                  tags$li("Compare microbial diversity across treatment groups"),
                  tags$li("Provide a reproducible, interactive analysis platform for stakeholders")
                )
              )
            )
          )
        )

      ),
      
      tabItem(tabName = "dashboard",

        # --- Dashboard header ---
        fluidRow(
          column(12,
            tags$div(class = "dash-header",
              tags$h3(icon("bacteria"), " chickMicro — Analysis Dashboard"),
              tags$p("Impact of Dulse supplementation on the poultry gut microbiome | AMR \u00b7 VF \u00b7 MGE profiling")
            )
          )
        ),

        # --- Single flat tabBox (no double nesting) ---
        fluidRow(
          tabBox(
            id = "dashboardTabs",
            width = 12,

            # --- Description tab ---
            tabPanel(
              title = tagList(icon("align-left"), " Overview"),
              fluidRow(
                column(6,
                  br(),
                  tags$img(src = "AMR_Pathogens.png", alt = "AMR Pathogens",
                           width = "100%", height = "auto",
                           style = "border-radius: 6px; box-shadow: 0 2px 8px rgba(0,0,0,.15);")
                ),
                column(6,
                  br(),
                  box(
                    width = 12, solidHeader = TRUE, status = "primary",
                    title = tagList(icon("leaf"), " Impact of Dulse on the Poultry Gut Microbiome"),
                    tags$p("The global spread of antimicrobial-resistant bacteria is a pressing public health crisis.
                            This study investigates ", tags$strong("Dulse"), " (", tags$em("Palmaria palmata"),
                           "), a protein-rich red seaweed, as a natural antibiotic alternative in broiler chicken production."),
                    tags$hr(style = "margin: 10px 0;"),
                    tags$ul(
                      style = "line-height: 1.9; padding-left: 18px;",
                      tags$li("WHO 2017 Priority Pathogen List: bacteria of highest urgency for new antibiotic development."),
                      tags$li("In 2019, ~1.3 million deaths (including 140,000 newborns) were directly attributed to AMR."),
                      tags$li("AMR-related deaths are projected to reach ", tags$strong("10 million/year by 2050"), "."),
                      tags$li("Poultry production is one of the largest agricultural users of antimicrobials worldwide.")
                    ),
                    tags$hr(style = "margin: 10px 0;"),
                    tags$p("This platform provides interactive tools to explore taxonomy, AMR genes, virulence factors
                            (VFs), and mobile genetic elements (MGEs) from metagenomic data generated across treatment groups.")
                  )
                )
              )
            ),

            # --- Data tab ---
            tabPanel(
              title = tagList(icon("upload"), " Data"),
              br(),
              csvFileUI("MDR_hvKp_dataSet")
            ),

            # --- Metadata tab ---
            tabPanel(
              title = tagList(icon("tags"), " Metadata"),
              br(),
              csvFileUI("MetadataSet")
            ),

            # --- Plot tab ---
            tabPanel(
              title = tagList(icon("chart-bar"), " Plot"),
              br(),
              MDRdataPlotUI("MDR_hvKp_dataPlot")
            ),

            # --- Stats tab ---
            tabPanel(
              title = tagList(icon("table"), " Stats Analysis"),
              br(),
              dataStatsUI("hosMicro_dataStat")
            )
          )
        )

      ),
      
      
      # <--- Contact Section ----->
      tabItem(tabName = "contact",

        fluidRow(
          column(12,
            tags$div(class = "dash-header",
              tags$h3(icon("envelope"), " Contact"),
              tags$p("Get in touch with the research team")
            )
          )
        ),

        fluidRow(
          box(
            title = tagList(icon("user"), " Researcher"), width = 4,
            status = "primary", solidHeader = TRUE,
            tags$p(tags$strong("Julio C. Ortega Cambara")),
            tags$p(icon("envelope"), tags$a("32104617@student.uwl.ac.uk",
                   href = "mailto:32104617@student.uwl.ac.uk")),
            tags$p(icon("graduation-cap"), " PhD Candidate — Computational Bioinformatics"),
            tags$p(icon("university"), " School of Biomedical Sciences"),
            tags$p(icon("map-marker-alt"), " University of West London"),
            tags$p(style = "color: #777; font-size: 12px;", "St Mary's Rd, London W5 5RF")
          ),

          box(
            title = tagList(icon("flask"), " Research Group"), width = 4,
            status = "info", solidHeader = TRUE,
            tags$p("This work is conducted within the ", tags$strong("School of Biomedical Sciences"),
                   " at the University of West London."),
            tags$p("Research focuses on the intersection of metagenomics, antimicrobial resistance,
                   and sustainable poultry production."),
            tags$p(icon("globe"), tags$a(" uwl.ac.uk",
                   href = "https://www.uwl.ac.uk", target = "_blank"))
          ),

          box(
            title = tagList(icon("code"), " App & Data"), width = 4,
            status = "success", solidHeader = TRUE,
            tags$p("Built with ", tags$strong("R Shiny"), " and ",
                   tags$strong("shinydashboard"), "."),
            tags$p("For bug reports, data queries, or collaboration requests,
                   please reach out via email."),
            tags$p(icon("database"), " Metagenomic data generated from broiler chicken gut samples
                   across three dietary treatment groups.")
          )
        )

      ),

      # <--- About Section ----->
      tabItem(tabName = "about",

        fluidRow(
          column(12,
            tags$div(class = "dash-header",
              tags$h3(icon("info-circle"), " About chickMicro"),
              tags$p("A web-based platform for poultry gut microbiome analysis")
            )
          )
        ),

        fluidRow(
          box(
            title = tagList(icon("bullseye"), " Purpose"), width = 6,
            status = "primary", solidHeader = TRUE,
            tags$p("chickMicro is a web-based application designed to explore the ",
                   tags$strong("poultry gut microbiome"), " — a critical tool for understanding
                   the complex microbial ecosystems within broiler chicken production."),
            tags$p("By providing detailed insights into the types and behaviours of microorganisms
                   present in poultry samples, the platform aids in:"),
            tags$ul(
              style = "line-height: 1.9; padding-left: 18px;",
              tags$li("Identifying potential sources of antimicrobial resistance (AMR)"),
              tags$li("Tracking the spread of antibiotic-resistant bacteria across treatment groups"),
              tags$li("Profiling virulence factors (VFs) and mobile genetic elements (MGEs)"),
              tags$li("Supporting targeted intervention strategies to reduce AMR in agriculture")
            )
          ),

          box(
            title = tagList(icon("layer-group"), " Tech Stack"), width = 3,
            status = "info", solidHeader = TRUE,
            tags$ul(
              style = "line-height: 2; padding-left: 18px; list-style: none;",
              tags$li(icon("r-project"),  " R 4.x"),
              tags$li(icon("desktop"),    " Shiny + shinydashboard"),
              tags$li(icon("chart-bar"), " ggplot2 + plotly"),
              tags$li(icon("table"),      " DT + gt + gtsummary"),
              tags$li(icon("dna"),        " vegan · VennDiagram · pheatmap"),
              tags$li(icon("file-alt"),   " rmarkdown reports")
            )
          ),

          box(
            title = tagList(icon("code-branch"), " Version"), width = 3,
            status = "warning", solidHeader = TRUE,
            tags$p(tags$strong("chickMicro v1.0")),
            tags$p(icon("calendar-alt"), " 2025"),
            tags$p(icon("university"),   " University of West London"),
            tags$p(icon("microscope"),   " School of Biomedical Sciences"),
            tags$hr(),
            tags$p(style = "font-size: 12px; color: #777;",
                   "For academic and research use. Data and results should be
                    interpreted in the context of the study design.")
          )
        )

      )
    )
  )
)

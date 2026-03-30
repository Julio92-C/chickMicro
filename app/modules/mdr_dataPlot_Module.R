# MDR dataPlot Module
# UI split into MDRreportUI (pre-loaded, no filters) and MDRexplorerUI (user-driven filters).
# Both share the same module server: MDRdataPlotServer.


# ---------------------------------------------------------------------------
# REPORT UI — auto-renders from pre-loaded data, no filter required
# ---------------------------------------------------------------------------
MDRreportUI <- function(id) {
  ns <- NS(id)

  tagList(

    # --- Summary boxes (full dataset) ---
    fluidRow(
      valueBoxOutput(ns("rpt_samples"), width = 2),
      valueBoxOutput(ns("rpt_taxa"),    width = 2),
      valueBoxOutput(ns("rpt_args"),    width = 2),
      valueBoxOutput(ns("rpt_vfs"),     width = 2),
      valueBoxOutput(ns("rpt_mges"),    width = 2),
      valueBoxOutput(ns("rpt_drugs"),   width = 2)
    ),

    tags$hr(style = "margin: 4px 0 18px 0;"),

    # --- Sidebar + main content ---
    fluidRow(

      # Left sidebar: section selector
      column(3,
        box(
          width = 12, solidHeader = TRUE, status = "primary",
          title = tagList(icon("layer-group"), " Plot Sections"),
          selectInput(ns("rpt_section"), label = NULL,
            choices = c(
              "Taxa Overview"        = "taxa",
              "Genetic Elements"     = "genes",
              "Microbiome Diversity" = "diversity",
              "Key Pathogens / AMR"  = "pathogens"
            ),
            width = "100%"
          ),
          tags$hr(),
          uiOutput(ns("rpt_desc"))
        )
      ),

      # Right main panel
      column(9,

        # --- Taxa Overview ---
        conditionalPanel(
          condition = sprintf("input['%s'] == 'taxa'", ns("rpt_section")),
          fluidRow(
            box(
              width = 12, solidHeader = TRUE, status = "teal",
              title = tagList(icon("bacterium"), " Taxa Overview"),
              fluidRow(
                column(7,
                  tags$p(tags$strong("Taxa count distribution by treatment group"),
                         style = "text-align:center; color:#555; margin-bottom:6px;"),
                  withSpinner(plotOutput(ns("taxaCount_violinPlot"), height = 420),
                              type = 6, color = "#1a73b8")
                ),
                column(5,
                  tags$p(tags$strong("Shared taxa across treatment groups"),
                         style = "text-align:center; color:#555; margin-bottom:6px;"),
                  withSpinner(plotOutput(ns("taxaCount_VennDiagram"), height = 420),
                              type = 6, color = "#1a73b8")
                )
              )
            )
          )
        ),

        # --- Genetic Elements ---
        conditionalPanel(
          condition = sprintf("input['%s'] == 'genes'", ns("rpt_section")),
          fluidRow(
            box(
              width = 12, solidHeader = TRUE, status = "purple",
              title = tagList(icon("dna"), " Genetic Elements Overview"),
              fluidRow(
                column(7,
                  tags$p(tags$strong("Gene count distribution by treatment group"),
                         style = "text-align:center; color:#555; margin-bottom:6px;"),
                  withSpinner(plotOutput(ns("ARGsCount_violinPlot"), height = 420),
                              type = 6, color = "#1a73b8")
                ),
                column(5,
                  tags$p(tags$strong("ARG / VF / MGE co-occurrence"),
                         style = "text-align:center; color:#555; margin-bottom:6px;"),
                  withSpinner(plotOutput(ns("geneticElements_occurrence"), height = 420),
                              type = 6, color = "#1a73b8")
                )
              )
            )
          )
        ),

        # --- Microbiome Diversity ---
        conditionalPanel(
          condition = sprintf("input['%s'] == 'diversity'", ns("rpt_section")),
          fluidRow(
            box(
              width = 12, solidHeader = TRUE, status = "navy",
              title = tagList(icon("circle-nodes"), " PCoA — Genetic Profile Ordination"),
              withSpinner(plotlyOutput(ns("PCoA_Analysis"), height = 420),
                          type = 6, color = "#1a73b8")
            )
          )
        ),

        # --- Key Pathogens / AMR ---
        conditionalPanel(
          condition = sprintf("input['%s'] == 'pathogens'", ns("rpt_section")),
          fluidRow(
            box(
              width = 12, solidHeader = TRUE, status = "danger",
              title = tagList(icon("bacteria"), " Key Pathogens — AMR Chord Diagram"),
              div(style = "padding: 15px;",
                withSpinner(plotOutput(ns("Key_Pathoges_AMR"), height = 800),
                            type = 6, color = "#1a73b8")
              ),
              tags$p(
                style = "font-size: 12px; color: #777; text-align: center; margin-top: 6px;",
                "Clockwise relationship: Samples \u2192 Taxa \u2192 AMR genes \u2192 Drug classes (Reference diet group)."
              )
            )
          )
        )

      ) # end column(9)
    )   # end fluidRow
  )     # end tagList
}


# ---------------------------------------------------------------------------
# EXPLORER UI — filter-driven, renders on Plot button click
# ---------------------------------------------------------------------------
MDRexplorerUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(

      # Left: filter panel
      column(3,
        box(
          width = 12, solidHeader = TRUE, status = "black", background = "gray",
          collapsible = TRUE,
          title = tagList(icon("filter"), " Filters"),
          selectInput(ns("DATABASE"), "Database:",
            choices  = c("card", "vfdb", "plasmidfinder"),
            selected = "card",
            multiple = TRUE
          ),
          tags$h5("Threshold filters:"),
          sliderInput(ns("NAME"),         "Min. taxa count", value = 20, min = 0, max = 100),
          sliderInput(ns("GENE"),         "Min. gene count", value = 10, min = 0, max = 100),
          sliderInput(ns("COVERAGE_PCT"), "Coverage (%)",    value = 90, min = 80, max = 100),
          sliderInput(ns("IDENTITY_PCT"), "Identity (%)",    value = 90, min = 80, max = 100),
          tags$hr(),
          fluidRow(
            column(6,
              actionButton(ns("plot_button"), tagList(icon("chart-bar"), " Plot"),
                style = "width:100%; background-color:#428bca; color:white; border:none;
                         border-radius:5px; padding:10px; font-size:14px; font-weight:bold;"
              )
            ),
            column(6,
              downloadButton(ns("download"), tagList(icon("file-alt"), " Report"),
                style = "width:100%; background-color:#428bca; color:white; border:none;
                         border-radius:5px; padding:10px; font-size:14px; font-weight:bold;"
              )
            )
          )
        )
      ),

      # Right: outputs
      column(9,
        fluidRow(
          infoBoxOutput(ns("general_stats"), width = 8),
          valueBoxOutput(ns("Drugs_count"),  width = 4)
        ),
        fluidRow(
          valueBoxOutput(ns("ARG_count"),  width = 4),
          valueBoxOutput(ns("VFs_count"),  width = 4),
          valueBoxOutput(ns("MGEs_count"), width = 4)
        ),
        fluidRow(
          tabBox(width = 12,

            tabPanel(
              title = tagList(icon("chart-bar"), " Frequency"),
              fluidRow(
                column(4),
                column(4,
                  selectInput(ns("mainVar_freq"), "Variable:",
                    choices = c("NAME", "GENE"), selected = "NAME", multiple = FALSE)
                ),
                column(4)
              ),
              withSpinner(plotlyOutput(ns("freq_plot"), height = 450),
                          type = 6, color = "#1a73b8")
            ),

            tabPanel(
              title = tagList(icon("chart-area"), " Relative Abundance"),
              withSpinner(plotlyOutput(ns("Relative_abundance"), height = 450),
                          type = 6, color = "#1a73b8")
            ),

            tabPanel(
              title = tagList(icon("chart-pie"), " Genetic Elements Profile"),
              fluidRow(
                column(6,
                  box(width = 12, solidHeader = TRUE, status = "warning",
                    title = tagList(icon("circle"), " Database distribution"),
                    plotOutput(ns("ARG_profile"), height = 420)
                  )
                ),
                column(6,
                  box(width = 12, solidHeader = TRUE, status = "warning",
                    title = tagList(icon("chart-bar"), " Elements by sample"),
                    plotOutput(ns("GeneticsElements_SAMPLE"), height = 420)
                  )
                )
              )
            ),

            tabPanel(
              title = tagList(icon("table-cells"), " Gut Microbiome"),
              fluidRow(
                column(4,
                  box(
                    width = 12, solidHeader = TRUE, status = "black", background = "gray",
                    title = tagList(icon("sliders"), " Heatmap Controls"),
                    sliderInput(ns("row_count"), "Row range:",
                                min = 1, max = 226, value = c(1, 50))
                  )
                ),
                column(8,
                  tags$p(style = "font-size: 12px; color: #777; padding: 6px 0;",
                    icon("info-circle"),
                    " Heatmap reflects the current filter selection. Click Plot to update."
                  )
                )
              ),
              fluidRow(
                withSpinner(plotOutput(ns("Gut_microbiome"), height = 750),
                            type = 6, color = "#1a73b8")
              )
            )

          )
        )
      )

    )
  )
}


# ---------------------------------------------------------------------------
# MODULE SERVER — shared by both Report and Explorer UIs
# ---------------------------------------------------------------------------
MDRdataPlotServer <- function(id, dataframe, metadata){
  stopifnot(is.reactive(dataframe))
  stopifnot(is.reactive(metadata))

  moduleServer(id, function(input, output, session){

    # -----------------------------------------------------------------------
    # REPORT: summary value boxes (full unfiltered dataset)
    # -----------------------------------------------------------------------
    output$rpt_samples <- renderValueBox({
      valueBox(
        value    = length(unique(dataframe()$SAMPLE)),
        subtitle = "Total Samples",
        icon     = icon("vial"),
        color    = "blue"
      )
    })

    output$rpt_taxa <- renderValueBox({
      valueBox(
        value    = length(unique(dataframe()$NAME)),
        subtitle = "Unique Taxa",
        icon     = icon("bacteria"),
        color    = "teal"
      )
    })

    output$rpt_args <- renderValueBox({
      n <- dataframe() %>% filter(DATABASE == "card") %>% distinct(GENE) %>% nrow()
      valueBox(value = n, subtitle = "AMR Genes",
               icon = icon("dna"), color = "purple")
    })

    output$rpt_vfs <- renderValueBox({
      n <- dataframe() %>% filter(DATABASE == "vfdb") %>% distinct(GENE) %>% nrow()
      valueBox(value = n, subtitle = "Virulence Factors",
               icon = icon("shield-virus"), color = "red")
    })

    output$rpt_mges <- renderValueBox({
      n <- dataframe() %>% filter(DATABASE == "plasmidfinder") %>% distinct(GENE) %>% nrow()
      valueBox(value = n, subtitle = "Mobile Genetic Elements",
               icon = icon("circle"), color = "orange")
    })

    output$rpt_drugs <- renderValueBox({
      n <- dataframe() %>%
        filter(DATABASE == "card") %>%
        separate_longer_delim(RESISTANCE, delim = ";") %>%
        distinct(RESISTANCE) %>%
        nrow()
      valueBox(value = n, subtitle = "Drug Classes",
               icon = icon("capsules"), color = "yellow")
    })

    output$rpt_desc <- renderUI({
      req(input$rpt_section)
      desc <- switch(input$rpt_section,
        "taxa"      = "Distribution and overlap of taxa across the three dietary treatment groups.",
        "genes"     = "Gene count by treatment and co-occurrence of ARGs, VFs, and MGEs.",
        "diversity" = "PCoA ordination of genetic profiles and gut microbiome presence/absence heatmap.",
        "pathogens" = "Chord diagram linking samples to key pathogens, AMR genes, and drug resistance classes."
      )
      tags$p(style = "font-size: 12px; color: #aaa; line-height: 1.6;", desc)
    })


    # -----------------------------------------------------------------------
    # EXPLORER: filter reactive (triggered by Plot button)
    # -----------------------------------------------------------------------
    filter_data <- eventReactive(input$plot_button, {
      validate(
        need(dataframe(), "Please input a data-sets as a csv file"),
        need(input$DATABASE, "Please select a database")
      )

      NAME_counts <- dataframe() %>%
        group_by(NAME) %>%
        summarize(count = n())

      gene_counts <- dataframe() %>%
        group_by(GENE) %>%
        summarize(count = n())

      dataframe() %>%
        filter(NAME %in% NAME_counts$NAME[NAME_counts$count > input$NAME],
               GENE %in% gene_counts$GENE[gene_counts$count > input$GENE],
               COVERAGE_PCT >= input$COVERAGE_PCT,
               IDENTITY_PCT >= input$IDENTITY_PCT,
               DATABASE %in% input$DATABASE)
    })

    # Download report
    output$download <- downloadHandler(
      filename = function() {
        paste("report-", Sys.Date(), ".html", sep = "")
      },
      content = function(file) {
        params <- list(filter_data   = filter_data(),
                       dataframe     = dataframe(),
                       metadata      = metadata(),
                       mainVar_freq  = input$mainVar_freq,
                       row_count     = input$row_count)

        rmarkdown::render("reports/shinyReport.Rmd", output_file = file,
                          params = params,
                          envir  = new.env(parent = globalenv()))
      }
    )

    # Download taxa violin PNG
    output$dl_taxaViolin <- downloadHandler(
      filename = function() paste0("taxa-violin-", Sys.Date(), ".png"),
      content = function(file) {
        abri_kraken2_merged <- merge(dataframe(), metadata(), by = "SAMPLE")
        abri_kraken2_filtered <- abri_kraken2_merged %>%
          arrange(TREATMENT) %>%
          group_by(TREATMENT, NAME) %>%
          summarise(Taxa_count = n(), .groups = "drop") %>%
          ungroup() %>%
          mutate(Taxa_log_count = log(Taxa_count + 1))
        anova_result <- aov(Taxa_log_count ~ TREATMENT, data = abri_kraken2_filtered)
        p <- ggplot(abri_kraken2_filtered, aes(x = TREATMENT, y = Taxa_log_count, fill = TREATMENT)) +
          geom_violin(trim = FALSE, scale = "width") +
          geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
          scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                                       "Seaweed"        = "#2ca02c",
                                       "Soyabean meal"  = "#ff7f0e")) +
          labs(x = "Treatment groups", y = "Log(Number of Taxa)") +
          theme_classic() +
          theme(legend.position = "none", text = element_text(size = 24)) +
          annotate("text", x = 2, y = max(abri_kraken2_filtered$Taxa_log_count) + 1.3,
                   label = paste("ANOVA p =", format(summary(anova_result)[[1]][["Pr(>F)"]][1], digits = 2)),
                   size = 5, color = "black")
        ggsave(file, plot = p, device = "png", width = 8, height = 6, dpi = 150)
      }
    )

    # Download Venn PNG
    output$dl_taxaVenn <- downloadHandler(
      filename = function() paste0("taxa-venn-", Sys.Date(), ".png"),
      content = function(file) {
        abri_kraken2_merged <- merge(dataframe(), metadata(), by = "SAMPLE")
        Wheat_soyabean_dataset <- abri_kraken2_merged %>%
          filter(TREATMENT == "Reference diet") %>% distinct() %>%
          rename(Wheat_soyabean = TREATMENT) %>% select(TAXID, Wheat_soyabean)
        Seaweed_dataset <- abri_kraken2_merged %>%
          filter(TREATMENT == "Seaweed") %>% distinct() %>%
          rename(Seaweed = TREATMENT) %>% select(TAXID, Seaweed)
        Soyabean_meal_dataset <- abri_kraken2_merged %>%
          filter(TREATMENT == "Soyabean meal") %>% distinct() %>%
          rename(Soyabean_meal = TREATMENT) %>% select(TAXID, Soyabean_meal)
        TREATMENT_dataset  <- merge(Wheat_soyabean_dataset, Seaweed_dataset, by = "TAXID", all = TRUE)
        TREATMENT_dataset1 <- merge(TREATMENT_dataset, Soyabean_meal_dataset, by = "TAXID", all = TRUE)
        TREATMENT_dataset1 <- TREATMENT_dataset1 %>%
          distinct(TAXID, .keep_all = TRUE) %>%
          remove_rownames() %>%
          column_to_rownames(var = "TAXID") %>%
          mutate(across(everything(), ~ if_else(. == "", NA_character_, .))) %>%
          filter(!if_all(everything(), is.na))
        TREATMENT_dataset_matrix        <- data.matrix(TREATMENT_dataset1, rownames.force = NA) %>% replace(is.na(.), 0)
        TREATMENT_dataset_matrix_Bmatrix <- as.matrix((TREATMENT_dataset_matrix > 0) + 0)
        sets <- apply(TREATMENT_dataset_matrix_Bmatrix, 2, function(col) which(col == 1))
        names(sets) <- colnames(TREATMENT_dataset_matrix_Bmatrix)
        venn.plot <- venn.diagram(
          x = sets,
          category.names = c("Wheat soyabean", "Seaweed", "Soyabean meal"),
          fill  = c(Wheat_soyabean = "#c2320e", Seaweed = "#04540a", Soyabean_meal = "#c2980e"),
          alpha = 0.3, height = 50, width = 50, filename = NULL, cat.cex = 1.2, cex = 2
        )
        png(file, width = 800, height = 800)
        grid.newpage(); grid.draw(venn.plot)
        dev.off()
      }
    )

    # -----------------------------------------------------------------------
    # EXPLORER: filter-dependent value boxes
    # -----------------------------------------------------------------------
    output$general_stats <- renderInfoBox({
      infoBox(
        title    = "General stats",
        value    = paste("Reads:", length(unique(filter_data()$SEQUENCE_ID)), "|",
                         "Samples:", length(unique(filter_data()$SAMPLE)), "|",
                         "Taxa:", length(unique(filter_data()$NAME)), "|",
                         "Genes:", length(unique(filter_data()$GENE)), "|",
                         "Resistance:", length(unique(filter_data()$RESISTANCE))),
        subtitle = "Kraken2 taxonomy | Abricate ARG/VF/MGE profiling",
        color    = "aqua",
        icon     = icon("chart-line"),
        fill     = TRUE
      )
    })

    output$ARG_count <- renderValueBox({
      arg_data <- filter_data() %>% filter(DATABASE == "card") %>% distinct(NAME, GENE)
      total    <- length(unique(arg_data$GENE))
      n_taxa   <- length(unique(filter_data()$NAME))
      pct      <- if (n_taxa > 0) round(100 * total / n_taxa, 2) else 0
      valueBox(value = paste0(total, " / ", pct, "%"),
               subtitle = "Antimicrobial Resistance Genes", icon("dna"), color = "purple")
    })

    output$VFs_count <- renderValueBox({
      vfs_data <- filter_data() %>% filter(DATABASE == "vfdb") %>% distinct(NAME, GENE)
      total    <- length(unique(vfs_data$GENE))
      n_taxa   <- length(unique(filter_data()$NAME))
      pct      <- if (n_taxa > 0) round(100 * total / n_taxa, 2) else 0
      valueBox(value = paste0(total, " / ", pct, "%"),
               subtitle = "Virulence Factors", icon("dna"), color = "purple")
    })

    output$MGEs_count <- renderValueBox({
      mges_data <- filter_data() %>% filter(DATABASE == "plasmidfinder") %>% distinct(NAME, GENE)
      total     <- length(unique(mges_data$GENE))
      n_taxa    <- length(unique(filter_data()$NAME))
      pct       <- if (n_taxa > 0) round(100 * total / n_taxa, 2) else 0
      valueBox(value = paste0(total, " / ", pct, "%"),
               subtitle = "Mobile Genetic Elements", icon("dna"), color = "purple")
    })

    output$Drugs_count <- renderValueBox({
      drugs_data <- filter_data() %>%
        filter(DATABASE == "card") %>%
        separate_longer_delim(RESISTANCE, delim = ";") %>%
        distinct(GENE, RESISTANCE)
      total  <- length(unique(drugs_data$RESISTANCE))
      n_taxa <- length(unique(filter_data()$NAME))
      pct    <- if (n_taxa > 0) round(100 * total / n_taxa, 2) else 0
      valueBox(value = paste0(total, " / ", pct, "%"),
               subtitle = "Drug Classes", icon("capsules"), color = "orange")
    })

    # Color palette
    palt1 <- as.character(paletteer_d("palettesForR::Cranes", n = 200))

    # -----------------------------------------------------------------------
    # EXPLORER: filter-dependent plots
    # -----------------------------------------------------------------------

    # Frequency bar chart
    output$freq_plot <- renderPlotly({
      a <- ggplot(filter_data(), aes(x = fct_infreq(.data[[input$mainVar_freq]]))) +
        geom_bar(fill = "skyblue") +
        geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
        theme_classic() +
        theme(axis.text   = element_text(size = 10.5),
              plot.title  = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 45)) +
        labs(x = "", y = "Count")
      plotly::ggplotly(a)
    })

    # Relative abundance stacked bar
    output$Relative_abundance <- renderPlotly({
      tryCatch({
        abri_kraken2_merged <- merge(filter_data(), metadata(), by = "SAMPLE")
        df_species <- abri_kraken2_merged %>%
          group_by(NAME) %>% arrange(NAME) %>%
          group_by(SAMPLE, NAME, TREATMENT) %>%
          summarise(Count = n(), .groups = "drop") %>%
          group_by(SAMPLE) %>%
          mutate(Percentage = Count / sum(Count) * 100)

        t          <- length(levels(as.factor(filter_data()$NAME)))
        t_palette  <- paletteer_d("palettesForR::Cranes", n = t)

        ra <- ggplot(df_species, aes(x = factor(SAMPLE), y = Percentage, fill = NAME)) +
          geom_bar(stat = "identity") +
          scale_fill_manual(values = as.character(t_palette)) +
          labs(x = "Sample groups", y = "Percentage", fill = "Taxonomic level") +
          theme_classic() +
          theme(legend.position   = "top",
                axis.text.x       = element_text(angle = 45, hjust = 1),
                text              = element_text(size = 14)) +
          facet_wrap(~ TREATMENT, scales = "free_x", nrow = 1)

        plotly::ggplotly(ra)
      }, error = function(e) {
        showNotification("Error generating plot. Number of colors exceeds palette limit (256).",
                         type = "error")
        NULL
      })
    })

    # ARG profile pie chart
    output$ARG_profile <- renderPlot({
      Gene_count <- table(filter_data()$DATABASE)
      pie(Gene_count,
          labels = paste0(names(Gene_count), " / ", Gene_count, " / ",
                          round(100 * Gene_count / sum(Gene_count), 2), "%"),
          col    = c("#b02d02", "#5c0404", "#c47a02"),
          cex    = 1.3, radius = 1, xlim = c(-1.5, 1.5))
    })

    # Genetic elements stacked bar by sample
    output$GeneticsElements_SAMPLE <- renderPlot({
      data_summary <- filter_data() %>%
        group_by(SAMPLE, DATABASE) %>%
        summarise(Count = n()) %>%
        group_by(SAMPLE) %>%
        mutate(Percentage = Count / sum(Count) * 100)

      ggplot(data_summary, aes(x = factor(SAMPLE), y = Percentage, fill = DATABASE)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("#b02d02", "#5c0404", "#c47a02")) +
        labs(x = "Sample groups by Treatment", y = "Percentage", fill = "Database") +
        theme_classic() +
        theme(legend.position   = "top",
              text              = element_text(size = 20),
              axis.text.x       = element_text(angle = 45, hjust = 1, size = 14),
              axis.text.y       = element_text(size = 14))
    })

    # -----------------------------------------------------------------------
    # REPORT: raw-data plots
    # -----------------------------------------------------------------------

    # Taxa Venn diagram by treatment
    output$taxaCount_VennDiagram <- renderPlot({
      abri_kraken2_merged <- merge(dataframe(), metadata(), by = "SAMPLE")

      Wheat_soyabean_dataset <- abri_kraken2_merged %>%
        filter(TREATMENT == "Reference diet") %>% distinct() %>%
        rename(Wheat_soyabean = TREATMENT) %>% select(TAXID, Wheat_soyabean)
      Seaweed_dataset <- abri_kraken2_merged %>%
        filter(TREATMENT == "Seaweed") %>% distinct() %>%
        rename(Seaweed = TREATMENT) %>% select(TAXID, Seaweed)
      Soyabean_meal_dataset <- abri_kraken2_merged %>%
        filter(TREATMENT == "Soyabean meal") %>% distinct() %>%
        rename(Soyabean_meal = TREATMENT) %>% select(TAXID, Soyabean_meal)

      TREATMENT_dataset  <- merge(Wheat_soyabean_dataset, Seaweed_dataset, by = "TAXID", all = TRUE)
      TREATMENT_dataset1 <- merge(TREATMENT_dataset, Soyabean_meal_dataset, by = "TAXID", all = TRUE)
      TREATMENT_dataset1 <- TREATMENT_dataset1 %>%
        distinct(TAXID, .keep_all = TRUE) %>%
        remove_rownames() %>% column_to_rownames(var = "TAXID") %>%
        mutate(across(everything(), ~ if_else(. == "", NA_character_, .))) %>%
        filter(!if_all(everything(), is.na))

      TREATMENT_dataset_matrix         <- data.matrix(TREATMENT_dataset1, rownames.force = NA) %>% replace(is.na(.), 0)
      TREATMENT_dataset_matrix_Bmatrix <- as.matrix((TREATMENT_dataset_matrix > 0) + 0)
      sets <- apply(TREATMENT_dataset_matrix_Bmatrix, 2, function(col) which(col == 1))
      names(sets) <- colnames(TREATMENT_dataset_matrix_Bmatrix)

      venn.plot <- venn.diagram(
        x = sets,
        category.names = c("Wheat soyabean", "Seaweed", "Soyabean meal"),
        fill  = c(Wheat_soyabean = "#c2320e", Seaweed = "#04540a", Soyabean_meal = "#c2980e"),
        alpha = 0.3, height = 50, width = 50, filename = NULL, cat.cex = 1.2, cex = 2
      )
      grid.newpage(); grid.draw(venn.plot)
    })

    # Taxa violin by treatment
    output$taxaCount_violinPlot <- renderPlot({
      tryCatch({
        abri_kraken2_merged <- merge(dataframe(), metadata(), by = "SAMPLE")
        abri_kraken2_filtered <- abri_kraken2_merged %>%
          arrange(TREATMENT) %>%
          group_by(TREATMENT, NAME) %>%
          summarise(Taxa_count = n(), .groups = "drop") %>%
          ungroup() %>%
          mutate(Taxa_log_count = log(Taxa_count + 1))

        anova_result <- aov(Taxa_log_count ~ TREATMENT, data = abri_kraken2_filtered)

        ggplot(abri_kraken2_filtered, aes(x = TREATMENT, y = Taxa_log_count, fill = TREATMENT)) +
          geom_violin(trim = FALSE, scale = "width") +
          geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
          scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                                       "Seaweed"        = "#2ca02c",
                                       "Soyabean meal"  = "#ff7f0e")) +
          labs(x = "Treatment groups", y = "Log(Number of Taxa)") +
          theme_classic() +
          theme(legend.position = "none", text = element_text(size = 16)) +
          annotate("text", x = 2, y = max(abri_kraken2_filtered$Taxa_log_count) + 1.3,
                   label = paste("ANOVA p =", format(summary(anova_result)[[1]][["Pr(>F)"]][1], digits = 2)),
                   size = 5, color = "black") +
          theme(text = element_text(size = 24))
      }, error = function(e) {
        showNotification("Error generating taxa violin plot. Adjust filter values.",
                         type = "error")
        NULL
      })
    })

    # Gene count violin by treatment
    output$ARGsCount_violinPlot <- renderPlot({
      tryCatch({
        abri_kraken2_merged <- merge(dataframe(), metadata(), by = "SAMPLE")
        abri_kraken2_filtered <- abri_kraken2_merged %>%
          arrange(NAME) %>%
          group_by(TREATMENT, GENE) %>%
          summarise(Gene_count = n(), .groups = "drop") %>%
          ungroup() %>%
          mutate(Gene_log_count = log(Gene_count + 1))

        anova_result <- aov(Gene_log_count ~ TREATMENT, data = abri_kraken2_filtered)

        ggplot(abri_kraken2_filtered, aes(x = TREATMENT, y = Gene_log_count, fill = TREATMENT)) +
          geom_violin(trim = FALSE, scale = "width") +
          geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
          scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                                       "Seaweed"        = "#2ca02c",
                                       "Soyabean meal"  = "#ff7f0e")) +
          labs(x = "Treatment groups", y = "Log(Number of Genes)") +
          theme_classic() +
          theme(legend.position = "none", text = element_text(size = 16),
                axis.text.y     = element_text(angle = 90, hjust = 0.5)) +
          annotate("text", x = 2, y = max(abri_kraken2_filtered$Gene_log_count) + 0.3,
                   label = paste("ANOVA p =", format(summary(anova_result)[[1]][["Pr(>F)"]][1], digits = 2)),
                   size = 5, color = "black") +
          theme(text = element_text(size = 24))
      }, error = function(e) {
        showNotification("Error generating gene violin plot. Adjust filter values.",
                         type = "error")
        NULL
      })
    })

    # ARG/VFs/MGEs co-occurrence Venn
    output$geneticElements_occurrence <- renderPlot({
      df_wide_ann_rows <- dataframe() %>%
        arrange(SAMPLE) %>%
        pivot_wider(names_from = DATABASE, values_from = GENE, values_fn = list) %>%
        group_by(NAME) %>%
        summarise(across(any_of(c("card", "vfdb", "plasmidfinder")),
                         ~ paste(na.omit(.), collapse = ", "))) %>%
        mutate(across(any_of(c("card", "vfdb", "plasmidfinder")),
                      ~ str_replace_all(., "(NuLL,|,NULL|,NULL,|NULL|, )", ""))) %>%
        na.omit() %>%
        remove_rownames() %>%
        column_to_rownames(var = "NAME") %>%
        mutate(across(everything(), ~ if_else(. == "", NA_character_, .)))

      SAMPLE2arg_data_matrix  <- data.matrix(df_wide_ann_rows, rownames.force = NA) %>% replace(is.na(.), 0)
      SAMPLE2arg_data_Bmatrix <- as.matrix((SAMPLE2arg_data_matrix > 0) + 0)
      sets <- apply(SAMPLE2arg_data_Bmatrix, 2, function(col) which(col == 1))
      names(sets) <- colnames(SAMPLE2arg_data_Bmatrix)

      venn.plot <- venn.diagram(
        x = sets,
        category.names = c("VFs", "AMR", "MGEs"),
        fill  = c(vfdb = "#EEAD0EFF", card = "#CD3333FF", plasmidfinder = "blue"),
        alpha = 0.5, height = 50, width = 50, filename = NULL, cat.cex = 1.5, cex = 2
      )
      grid.newpage(); grid.draw(venn.plot)
    })

    # PCoA
    output$PCoA_Analysis <- renderPlotly({
      suppressWarnings({
        abri_kraken2_merged <- merge(dataframe(), metadata(), by = "SAMPLE")
        abri_kraken2_arranged <- abri_kraken2_merged %>%
          arrange(GENE) %>%
          select(SAMPLE, NAME, GENE, TREATMENT)

        df_wide_ann_rows <- abri_kraken2_arranged %>%
          arrange(SAMPLE) %>%
          pivot_wider(names_from = SAMPLE, values_from = c(NAME), values_fn = length) %>%
          group_by(GENE) %>%
          summarise(across(everything(), ~ paste(., collapse = ", "))) %>%
          mutate(across(-GENE, ~ str_replace_all(., "(NA,|,NA|,NA,|NA|,| )", "")))

        df_wide_ann_rows <- df_wide_ann_rows %>%
          na.omit() %>%
          remove_rownames() %>%
          column_to_rownames(var = "GENE") %>%
          mutate(across(everything(), ~ if_else(. == "", NA_character_, .)))

        df_clean <- df_wide_ann_rows[!apply(df_wide_ann_rows, 1, function(row) all(row == "" | is.na(row))), ]

        SAMPLE2arg_data_matrix  <- data.matrix(df_clean, rownames.force = NA) %>% replace(is.na(.), 0)
        bray_curtis_dist        <- vegdist(SAMPLE2arg_data_matrix, method = "bray")
        pcoa_result             <- cmdscale(bray_curtis_dist, eig = TRUE, k = 2)
        eig_values              <- pcoa_result$eig
        var_explained           <- eig_values / sum(eig_values) * 100

        pcoa_df <- data.frame(GENE = rownames(SAMPLE2arg_data_matrix),
                              PC1  = pcoa_result$points[, 1],
                              PC2  = pcoa_result$points[, 2])

        TREATMENTData <- abri_kraken2_arranged %>%
          select(SAMPLE, GENE, TREATMENT) %>% distinct()

        ann_pcoa       <- inner_join(pcoa_df, TREATMENTData, by = "GENE")
        TREATMENT_col  <- as.factor(ann_pcoa$TREATMENT)
        l              <- length(levels(TREATMENT_col))
        l_palette      <- paletteer_d("ggsci::default_igv", n = l)

        pcaoa_plot <- ggplot(ann_pcoa, aes(x = PC1, y = PC2, colour = TREATMENT)) +
          geom_point(size = 3) +
          stat_ellipse(type = "norm", lwd = 0.8) +
          scale_color_manual(values = l_palette) +
          labs(color = "Treatment",
               x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
               y = paste0("PC2 (", round(var_explained[2], 2), "%)")) +
          theme_classic() +
          theme(text = element_text(size = 16)) +
          scale_x_continuous(limits = c(-0.6, NA)) +
          scale_y_continuous(limits = c(-0.8, NA))

        plotly::ggplotly(pcaoa_plot)
      })
    })

    # Key Pathogens chord diagram
    output$Key_Pathoges_AMR <- renderPlot({
      abri_kraken2_merged <- merge(dataframe(), metadata(), by = "SAMPLE")

      modify_strings <- function(df, column_NAME, target_strings, replacement_strings) {
        for (i in seq_along(target_strings)) {
          df <- df %>%
            mutate(!!sym(column_NAME) := str_replace_all(!!sym(column_NAME),
                                                          target_strings[i], replacement_strings[i]))
        }
        return(df)
      }

      target_strings      <- c("carbapenem;cephalosporin;penam")
      replacement_strings <- c("Beta-lactam")
      abri_kraken2_filtered <- modify_strings(abri_kraken2_merged, "RESISTANCE",
                                              target_strings, replacement_strings)

      abri_kraken2_filtered <- abri_kraken2_filtered %>%
        filter(grepl("card", DATABASE)) %>%
        mutate(GENE = case_when(
          grepl("vanW_gene_in_vanB_cluster",              GENE) ~ "vanWB",
          grepl("vanR_gene_in_vanB_cluster",              GENE) ~ "vanRB",
          grepl("vanX_gene_in_vanB_cluster",              GENE) ~ "vanXB",
          grepl("vanY_gene_in_vanB_cluster",              GENE) ~ "vanYB",
          grepl("vanS_gene_in_vanB_cluster",              GENE) ~ "vanSB",
          grepl("vanH_gene_in_vanB_cluster",              GENE) ~ "vanHB",
          grepl("Escherichia_coli_ampC_beta-lactamase",   GENE) ~ "ampC",
          grepl("Escherichia_coli_mdfA",                  GENE) ~ "mdfA",
          grepl("Escherichia_coli_emrE",                  GENE) ~ "emrE",
          TRUE ~ GENE
        ),
        GENE       = ifelse(grepl("AAC(6')-Ie-APH(2'')-Ia_bifunctional_protein", GENE),
                            "aac(6')-Ie/aph(2'')-Ia", GENE),
        NAME       = ifelse(grepl("Chryseobacterium panacisoli", NAME), "C. panacisoli", NAME),
        RESISTANCE = ifelse(nchar(RESISTANCE) > 20, "Multi-drug", RESISTANCE)) %>%
        mutate(RESISTANCE = str_replace_all(RESISTANCE, "_", " ") %>% str_to_title()) %>%
        arrange(NAME)

      pathogen_AMR <- abri_kraken2_filtered %>%
        filter(!grepl("Alistipes|Azorhizobium|cellular organisms|Bacillota|Clostridia|Lachnospiraceae|Bacteria|root|Enterobacterales|Enterobacteriaceae|Terrabacteria group|Homo sapiens|Pseudomonadota|Bacteroidales|Bifidobacterium|Bacteroidota", NAME)) %>%
        filter(!(str_count(NAME, "\\S+") == 1 & NAME != "Enterococcus")) %>%
        filter(TREATMENT == "Reference diet") %>%
        select(SAMPLE, NAME, GENE, RESISTANCE) %>%
        mutate(NAME = str_replace(NAME, "^(\\w)\\w+\\s(\\w+)", "\\1. \\2")) %>%
        distinct()

      pathogen_AMR[64, 2] <- "C. innocuum"
      pathogen_AMR[65, 2] <- "C. scindens"
      pathogen_AMR[67:68, 2] <- "u. Subdoli."
      pathogen_AMR[66, 2] <- "u. Plantact."

      df         <- pathogen_AMR
      categories <- unique(c(sort(df$SAMPLE), sort(df$NAME), sort(df$GENE), sort(df$RESISTANCE)))
      mat        <- matrix(0, nrow = length(categories), ncol = length(categories),
                           dimnames = list(categories, categories))
      for (i in 1:nrow(df)) {
        mat[df$SAMPLE[i], df$NAME[i]]       <- 1
        mat[df$NAME[i],   df$GENE[i]]       <- 1
        mat[df$GENE[i],   df$RESISTANCE[i]] <- 1
      }

      color_sectors <- c(
        paletteer_d("ggsci::default_igv", n = length(unique(df$SAMPLE))),
        paletteer_d("ggsci::default_igv", n = length(unique(df$NAME))),
        paletteer_d("ggsci::default_igv", n = length(unique(df$GENE))),
        paletteer_d("ggsci::default_igv", n = length(unique(df$RESISTANCE)))
      )

      circos.clear()
      circos.par(track.height = 0.1, start.degree = 115, gap.degree = 2,
                 canvas.xlim = c(-1, 1), canvas.ylim = c(-1, 1),
                 circle.margin = c(1, 1), unit.circle.segments = 500)
      chordDiagram(mat, transparency = 0.5, annotationTrack = "grid",
                   scale = FALSE, directional = 1, diffHeight = mm_h(3),
                   grid.col = color_sectors,
                   preAllocateTracks = list(track.height = 0.1,
                                            unit.circle.segments = 500,
                                            start.degree = 95, scale = TRUE))
      circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$ylim[1],
                    get.cell.meta.data("sector.index"),
                    facing = "clockwise", niceFacing = TRUE,
                    adj = c(0, 0.8), cex = 1)
      }, bg.border = NA)
    })

    # Gut microbiome heatmap (reactive — uses filter_data)
    output$Gut_microbiome <- renderPlot({
      validate(
        need(!is.null(filter_data()), "Run a filter first — set your thresholds and click Plot."),
        need(nrow(filter_data()) > 0,  "No data after current filters.")
      )
      abri_kraken2_merged <- merge(filter_data(), metadata(), by = "SAMPLE")
      abri_kraken2_clean  <- abri_kraken2_merged %>%
        filter(!(str_count(NAME, "\\S+") == 1 & NAME != "Enterococcus"))

      df_wide <- abri_kraken2_clean %>%
        arrange(SAMPLE) %>%
        select(SAMPLE, TAXID, NAME) %>%
        pivot_wider(names_from = SAMPLE, values_from = TAXID, values_fn = length) %>%
        group_by(NAME) %>%
        summarise(across(everything(), ~ paste(., collapse = ", "))) %>%
        mutate(across(-NAME, ~ str_replace_all(., "(NA,|,NA|,NA,|NA|,| )", "")))

      processed_mgs2arg_data <- df_wide %>%
        remove_rownames() %>%
        column_to_rownames(var = "NAME") %>%
        mutate(across(everything(), ~ na_if(., "")))

      SAMPLE2arg_data <- processed_mgs2arg_data %>%
        filter(rowSums(is.na(.)) < ncol(.))

      SAMPLE2arg_data_matrix  <- data.matrix(SAMPLE2arg_data, rownames.force = NA) %>% replace(is.na(.), 0)
      SAMPLE2arg_data_Bmatrix <- as.matrix((SAMPLE2arg_data_matrix > 0) + 0)

      rearrange_matrix <- function(mat) mat[order(rowSums(mat), decreasing = TRUE), ]
      sorted_matrix    <- rearrange_matrix(SAMPLE2arg_data_Bmatrix)

      callback <- function(hc, mat) {
        sv   <- svd(t(mat))$v[, 1]
        dend <- reorder(as.dendrogram(hc), wts = sv)
        as.hclust(dend)
      }

      SAMPLE2arg_data_transposed <- as.data.frame(t(SAMPLE2arg_data)) %>%
        mutate(SAMPLE = rownames(.)) %>% select(SAMPLE)

      MetadataLocations <- abri_kraken2_clean %>%
        select(SAMPLE, TREATMENT) %>% distinct()

      ann_col <- inner_join(SAMPLE2arg_data_transposed, MetadataLocations, by = "SAMPLE") %>%
        na.omit() %>% remove_rownames() %>% column_to_rownames(var = "SAMPLE")

      d_palette  <- paletteer_d("ggsci::category10_d3", n = nlevels(as.factor(ann_col$TREATMENT)))
      df4        <- ann_col %>% distinct(TREATMENT) %>% mutate(color = d_palette)
      type_colors <- setNames(as.character(df4$color), df4$TREATMENT)

      df_wide_ann_rows <- abri_kraken2_clean %>%
        arrange(SAMPLE) %>%
        filter(grepl("card|vfdb|plasmidfinder", DATABASE)) %>%
        pivot_wider(names_from = DATABASE, values_from = GENE, values_fn = list) %>%
        group_by(NAME) %>%
        summarise(across(any_of(c("card", "vfdb", "plasmidfinder")),
                         ~ paste(., collapse = ", "))) %>%
        mutate(across(any_of(c("card", "vfdb", "plasmidfinder")),
                      ~ str_replace_all(., "(NULL,|,NULL|,NULL,|NULL| )", ""))) %>%
        na.omit() %>% remove_rownames() %>% column_to_rownames(var = "NAME") %>%
        mutate(across(everything(), ~ if_else(. == "", NA_character_, .)))

      classify <- function(input_string, type) {
        if (is.na(input_string)) return(NA)
        n <- length(unique(unlist(strsplit(input_string, ","))))
        switch(type,
          virulence = ifelse(n == 1, "Virulent",      "Hypervirulent"),
          amr       = ifelse(n == 1, "Drug-resistant","Multidrug-resistant"),
          MGE       = ifelse(n == 1, "Single MGE",    "Multi-MGEs")
        )
      }

      ann_rows <- df_wide_ann_rows %>%
        mutate(VFs  = sapply(vfdb,          classify, type = "virulence"),
               AMR  = sapply(card,          classify, type = "amr"),
               MGEs = sapply(plasmidfinder, classify, type = "MGE")) %>%
        select(VFs, AMR, MGEs)

      clean_string <- function(x) str_to_title(str_replace_all(x, "_", " "))

      Resistance_category <- abri_kraken2_clean %>%
        distinct(NAME, GENE, RESISTANCE) %>%
        mutate(RESISTANCE = sapply(RESISTANCE, clean_string)) %>%
        group_by(NAME) %>%
        summarise(RESISTANCE = paste(unique(RESISTANCE), collapse = ", ")) %>%
        column_to_rownames(var = "NAME") %>%
        mutate(across(1, ~ str_replace_all(., "(NA,|,NA|,NA,|NA|,| )", ""))) %>%
        mutate(RESISTANCE = ifelse(nchar(RESISTANCE) > 20, "Multi-drug", RESISTANCE)) %>%
        mutate(across(everything(), ~ if_else(. == "", NA_character_, .))) %>%
        select(RESISTANCE)

      ann_row_all <- merge(ann_rows, Resistance_category, by = 0) %>%
        rename(DRUG = RESISTANCE) %>%
        column_to_rownames(var = "Row.names") %>%
        select(DRUG, AMR, VFs, MGEs)

      r_palette        <- paletteer_d("ggsci::default_igv", n = length(unique(ann_row_all$DRUG)) - 1)
      ann_row_drug_col <- ann_row_all %>% select(DRUG) %>% na.omit() %>% distinct(DRUG) %>%
        mutate(color = r_palette)
      drug_colors <- setNames(as.character(ann_row_drug_col$color), ann_row_drug_col$DRUG)

      ann_colors <- list(
        TREATMENT = type_colors,
        DRUG      = drug_colors,
        AMR       = c(`Drug-resistant` = "#a02d02", `Multidrug-resistant` = "#522501"),
        MGEs      = c(`Multi-MGEs` = "#5c0404",     `Single MGE` = "#fcd2d2"),
        VFs       = c(`Virulent` = "#fce5d2",        `Hypervirulent` = "#c47a02")
      )

      tryCatch({
        pheatmap(sorted_matrix[input$row_count[1]:input$row_count[2], ],
                 display_numbers = FALSE, cluster_cols = TRUE, cluster_rows = FALSE,
                 scale = "none", clustering_callback = callback,
                 border_color = "NA", color = c("#CCCCCCFF", "#666666FF"),
                 legend_breaks = c(0, 1), legend_labels = c("Absent", "Present"),
                 annotation_row = ann_row_all, annotation_col = ann_col[1],
                 show_rownames = TRUE, annotation_colors = ann_colors,
                 fontsize_row = 14, fontsize_col = 14, fontsize = 14)
      }, error = function(e) {
        message("Heatmap error: ", e$message)
      })
    })

    # Return filtered data for Stats module
    return(filter_data)
  })
}

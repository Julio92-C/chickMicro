# MDR-hvKp dataPlot Module


# Module UI function
MDRdataPlotUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(3,
             box(width = 12, 
                 solidHeader = TRUE, 
                 collapsible=TRUE,
                 status = "black",
                 title = "Choose a Database:",
                 background = "gray",
                 
                 selectInput(ns("DATABASE"),
                             label = "Database",
                             choices = c("card",
                                         "vfdb",
                                         "plasmidfinder"
                             ),
                             selected = "card",
                             multiple = TRUE
                 ),
                 
                 h4("Filter data-set by:"),
                 sliderInput(ns("NAME"),
                             label = "Taxa count",
                             value = 20, min = 0, max = 100
                             
                 ),
                 sliderInput(ns("GENE"),
                             label = "Gene count",
                             value = 10, min = 0, max = 100
                             
                 ),
                 
                 sliderInput(ns("COVERAGE_PCT"),
                             label = "Gene Coverage (%)",
                             value= 90, min = 80, max = 100

                             ),
                 sliderInput(ns("IDENTITY_PCT"),
                             label = "Gene Identity (%)",
                             value = 90, min = 80, max = 100

                 ),
                 
                 
                 
                 fluidRow(
                   column(6,
                          actionButton(ns("plot_button"), icon = icon("chart-bar"), "Plot", 
                                       style = "
                                       width: 100%;
                                       background-color: #428bca;
                                       color: white;
                                       border: none;
                                       border-radius: 5px;
                                       padding: 10px 20px;
                                       font-size: 16px;
                                       font-weight: bold;"
                          )
                   ),
                   column(6,
                          downloadButton(ns("download"), "Report",
                                         style = "
                                         width: 100%;
                                         background-color: #428bca;
                                         color: white;
                                         border: none;
                                         border-radius: 5px;
                                         padding: 10px 20px;
                                         font-size: 16px;
                                         font-weight: bold;"
                          )  
                   )
                 )
                 
                 
                 
                 
             )
      ),
      
      column(9,
             fluidRow(
               
               infoBoxOutput(ns("general_stats"), width = 8),
               
               valueBoxOutput(ns("Drugs_count"), width = 4)
               
               
             ),
             
             fluidRow(
               valueBoxOutput(ns("ARG_count"), width = 4),
               
               valueBoxOutput(ns("VFs_count"), width = 4),
               
               valueBoxOutput(ns("MGEs_count"), width = 4)
               
               
               
               
             ),
             
             fluidRow(
               
               tabBox(width = 12,
                      height = 1200,
                      
                      tabPanel("Taxa Distribution",
                               fluidRow(
                                 
                                   plotlyOutput(ns("freq_plot"), height = 500)
                                 
                               ),
                               
                               fluidRow(
                                 column(4,
                                        
                                 ),
                                 column(4,
                                       selectInput(ns("mainVar_freq"),
                                                   label = "Choose a genotype:",
                                                   choices = c("NAME", "GENE"),
                                                   selected = "NAME",
                                                   multiple = FALSE,
                                       )
                                 ),
                                 
                                 column(4,
                                   
                                 )
                               ),
                               
                               fluidRow(
                                 column(2,
                                        
                                 ),
                                 column(8,
                                        plotOutput(ns("taxaCount_violinPlot"), height = 500),
                                 ),
                                 column(2,
                                        
                                 ),
                               ),
                               
                               
                               
                               
                      ),
                      
                      
                      tabPanel("Relative Abundance",
                               fluidRow(
                                 plotlyOutput(ns("Relative_abundance"), height = 500)
                               ),
                               
                               fluidRow(
                                  column(2,
                                         
                                  ),
                                  column(8,
                                         h4("")
                                  ),
                                  column(2,
                                         
                                  ),
                                 ),
                               
                               fluidRow(
                                 column(2,
                                        
                                 ),
                                 column(8,
                                        plotOutput(ns("taxaCount_VennDiagram"), height = 600, width = 600)
                                 ),
                                 column(2,
                                        
                                 ),
                                 
                                 
                               )
                               
                      ),
                      
                      
                      
                      tabPanel("ARG/VFs/MGEs profile",
                               fluidRow(
                                 plotOutput(ns("ARG_profile"), height = 500),
                               ),
                               
                               fluidRow(
                                 plotOutput(ns("GeneticsElements_SAMPLE"), height = 500)
                               )
                               
                               

                      ),
                      
                      tabPanel("ARG/VFs/MGEs co-distribution",
                               fluidRow(
                                 plotlyOutput(ns("PCoA_Analysis"), height = 500),
                               ),
                               
                               fluidRow(
                                 column(6,
                                        plotOutput(ns("ARGsCount_violinPlot"), height = 500)  
                                 ),
                                 column(6,
                                        plotOutput(ns("geneticElements_occurrence"), height = 500, width = 500)


                                 )
                                 
                               )
                                 
                               
                      ),
                       
                      tabPanel("Gut microbiome",
                               fluidRow(
                                 plotOutput(ns("Gut_microbiome"), height = 900),
                               ),
                               
                               fluidRow(
                                 column(3,
                                        
                                 ),
                                 column(6,
                                        box(width = 12, 
                                            solidHeader = TRUE, 
                                            collapsible=TRUE,
                                            status = "black",
                                            title = "Dynamic Row Control for Heatmap:",
                                            background = "gray",
                                            # Slider input for selecting a range of values
                                            sliderInput(ns("row_count"), 
                                                "Select a range:",
                                                 min = 1, 
                                                 max = 226, 
                                                 value = c(1, 50))
                                    
                                        )
                                        
                                        
                                 ),
                                 column(3,
                                        
                                 )
                               ),
                               
                               
                      ),
                      
                      tabPanel("Key Pathogens/AMR",
                            fluidRow(
                                 div(
                                   style = "padding: 30px;",  # Adjust padding as needed
                                   plotOutput(ns("Key_Pathoges_AMR"), height = "800px")
                                 )
                               ),
                               fluidRow(
                                 h4("Chord diagram highlighting the clockwise relationship between the nodes: sample, taxa, ARGs and drug classes.")
                               )
                               
                          )
                     )
                      
               )
               
             )
      )
      
    )
  
  
}

# Module server function
MDRdataPlotServer <- function(id, dataframe, metadata){
  stopifnot(is.reactive(dataframe))
  stopifnot(is.reactive(metadata))
  
  moduleServer(id, function(input, output, session){
    
    # Filter data by user input
    filter_data <- eventReactive(input$plot_button, {
      validate(
        need(dataframe(), "Please input a data-sets as a csv file"),
        need(input$DATABASE, "Please select a database")
      )
      
      # Group the data by NAME and count the occurrences
      NAME_counts <- dataframe() %>% 
        group_by(NAME) %>%
        summarize(count = n())
      
      gene_counts <- dataframe() %>% 
        group_by(GENE) %>%
        summarize(count = n())
      
      # print(gene_counts)


      filter_data <- dataframe() %>%
        filter(NAME %in% NAME_counts$NAME[NAME_counts$count > input$NAME],
               GENE %in% gene_counts$GENE[gene_counts$count > input$GENE],
               COVERAGE_PCT >= input$COVERAGE_PCT,
               IDENTITY_PCT >= input$IDENTITY_PCT,
               DATABASE %in% input$DATABASE)
      
    })
    
    # Download plots to a HTML file
    output$download <- downloadHandler(
      filename = function() {
        paste("report-", Sys.Date(), ".html", sep="")
      },


      content = function(file) {
        # Code to generate the report goes here
        
        # Set up parameters to pass to Rmd document
        params <- list(filter_data = filter_data(),
                       dataframe = dataframe(),
                       metadata = metadata(),
                       mainVar_freq = input$mainVar_freq,
                       row_count = input$row_count)
        

        rmarkdown::render("reports/shinyReport.Rmd", output_file = file, 
                          params = params, 
                          envir = new.env(parent = globalenv())
                          )
      }
    )
    
    # Display an infoBox with general stats
    output$general_stats <- renderInfoBox({
      infoBox(
        title = "General stats",
        value = paste("No. of Reads:", length(unique(filter_data()$SEQUENCE_ID)), "|",
                      "Samples:", length(unique(filter_data()$SAMPLE)), "|",
                      "Taxa:", length(unique(filter_data()$NAME)), "|",
                      "Genes:", length(unique(filter_data()$GENE)), "|",
                      "Resistance:", length(unique(filter_data()$RESISTANCE)) 
                      
                      ),
        subtitle = "Taxonomy classification carried out using Kraken2. Profiling of ARGs, VFs, and MGEs with Abricate",
        color = "aqua",
        icon = icon("chart-line"),
        fill = TRUE
      )
    })
    
    output$ARG_count <- renderValueBox({
      arg_data <- filter_data() %>%
        filter(DATABASE == "card") %>%
        distinct(NAME, GENE) %>%
        group_by(GENE) %>%
        arrange(GENE)
      
      total_ARG_count <- length(unique(arg_data$GENE))
     
      
      valueBox(
        value = paste0(total_ARG_count, "/",
                       round(100 * total_ARG_count/length(unique(filter_data()$NAME)), 2), "%"),
        subtitle = "Total Antimacrobial Resistance genes", icon("dna"),
        color = "purple"
      )
    })
    
    
    
    output$VFs_count <- renderValueBox({
      vfs_data <- filter_data() %>%
        filter(DATABASE == "vfdb") %>%
        distinct(NAME, GENE) %>%
        group_by(GENE) %>%
        arrange(GENE)
      
      total_VFs_count <- length(unique(vfs_data$GENE))
      
      
      valueBox(
        value = paste0(total_VFs_count, "/",
                       round(100 * total_VFs_count/length(unique(filter_data()$NAME)), 2), "%"),
        subtitle = "Total Virulence Factors", icon("dna"),
        color = "purple"
      )
    })
    
    output$MGEs_count <- renderValueBox({
      mges_data <- filter_data() %>%
        filter(DATABASE == "plasmidfinder") %>%
        distinct(NAME, GENE) %>%
        group_by(GENE) %>%
        arrange(GENE)
      
      total_mges_count <- length(unique(mges_data$GENE))
      
      
      valueBox(
        value = paste0(total_mges_count, "/",
                       round(100 * total_mges_count/length(unique(filter_data()$NAME)), 2), "%"),
        subtitle = "Total Mobile Genetic Elements", icon("dna"),
        color = "purple"
      )
    })
    
    output$Drugs_count <- renderValueBox({
      drugs_data <- filter_data() %>%
        filter(DATABASE == "card") %>%
        separate_longer_delim(RESISTANCE, delim = ";") %>%
        distinct(GENE, RESISTANCE) %>%
        group_by(GENE) %>%
        arrange(GENE) 
         
      # drugs_class <- drugs_data$RESISTANCE
      # print(drugs_class)
      
      total_drugs_count <- length(unique(drugs_data$RESISTANCE))
      
      valueBox(
        value = paste0(total_drugs_count, "/",
                       round(100 * total_drugs_count/length(unique(filter_data()$NAME)), 2), "%"),
        subtitle = "Total Drug classes", icon("capsules"),
        color = "orange"
      )
    })
    
    # Color palettes
    palt1 <- paletteer_d("palettesForR::Cranes", n = 200)
    palt1 <- as.character(palt1)
    
    
    
    # Render ST type frequency plot based on the filter data sets
    output$freq_plot <- renderPlotly({
      # Plot ST type frequency
      a <- ggplot(filter_data(), aes(x= fct_infreq(.data[[input$mainVar_freq]]),
                                     ))
      a <- a + geom_bar(fill = "skyblue")
      a <- a + theme_classic()
      a <- a + geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5)  # Add count labels on top of the bars
      a <- a + labs(fill = "Legend")
      a <- a + theme(legend.position = "Bottom")
      a <- a + theme(legend.justification = "center")
      a <- a + theme(axis.text = element_text(size = 10.5))
      a <- a + theme(plot.title = element_text(hjust = 0.5))
      a <- a + theme(axis.text.x = element_text(angle = 45))
      a <- a + labs(x="", y="Count")
      a
      
      # Plot Model type frequency with PLOTLY ####
      ST_frequency_plot <- plotly::ggplotly(a)
      ST_frequency_plot
      
    })
    
    
    # Render TREATMENT venn Diagram based on the taxa count
    output$taxaCount_VennDiagram <- renderPlot({
      # Merge abri_kraken2_filtered and metadata dataframe
      abri_kraken2_merged <- merge(dataframe(), metadata(), by ="SAMPLE")
      
      # Wheat-soyabean dataset
      Wheat_soyabean_dataset <- abri_kraken2_merged %>%
        filter(TREATMENT == "Reference diet") %>%
        distinct() %>%
        rename(Wheat_soyabean = TREATMENT) %>%
        select(TAXID, Wheat_soyabean)
      
      # Seaweed dataset
      Seaweed_dataset <- abri_kraken2_merged %>%
        filter(TREATMENT == "Seaweed") %>%
        distinct() %>%
        rename(Seaweed = TREATMENT) %>%
        select(TAXID, Seaweed)
      
      # Soyabean meal dataset
      Soyabean_meal_dataset <- abri_kraken2_merged %>%
        filter(TREATMENT == "Soyabean meal") %>%
        distinct() %>%
        rename(Soyabean_meal = TREATMENT) %>%
        select(TAXID, Soyabean_meal)
      
      # Merge all TREATMENTs datasets
      TREATMENT_dataset <- merge(Wheat_soyabean_dataset, Seaweed_dataset, by="TAXID", all = TRUE)
      TREATMENT_dataset1 <- merge(TREATMENT_dataset, Soyabean_meal_dataset, by="TAXID", all = TRUE)
      #colnames(TREATMENT_dataset1)
      
      # Assign the values of col1 as row NAMEs
      TREATMENT_dataset1 <- TREATMENT_dataset1 %>%
        distinct(TAXID, .keep_all = TRUE) %>%
        remove_rownames %>% 
        column_to_rownames(var="TAXID") %>%
        mutate(across(everything(), ~ if_else(. == "", NA_character_, .))) %>%
        filter(!if_all(everything(), is.na))
      
      # Convert from a dataframe to a matrix
      TREATMENT_dataset_matrix <- data.matrix(TREATMENT_dataset1, rownames.force = NA)
      
      # Count NA in data sets
      sum(is.na(TREATMENT_dataset_matrix))
      
      # Replace all NA values with 0 values
      TREATMENT_dataset_matrix <- TREATMENT_dataset_matrix %>% replace(is.na(.), 0)
      
      # Convert to a binary matrix
      TREATMENT_dataset_matrix_Bmatrix <- as.matrix((TREATMENT_dataset_matrix > 0) + 0)
      
      # Create a list of sets
      sets <- apply(TREATMENT_dataset_matrix_Bmatrix, 2, function(col) which(col == 1))
      names(sets) <- colnames(TREATMENT_dataset_matrix_Bmatrix)
      print(sets)
      class(sets)
      
      ### Plot the Venn diagram
      venn.plot <- venn.diagram(
        x = sets,
        category.names = c("Wheat soyabean", "Seaweed", "Soyabean meal"),
        fill = c(Wheat_soyabean="#c2320e", Seaweed="#04540a", Soyabean_meal="#c2980e"), 
        alpha = 0.3,
        height = 50, width = 50,
        filename = NULL,
        cat.cex = 1.2,  # Adjust this value to change the font size
        cex = 2       # Adjust this value to change the font size of the numbers
      )
      grid.newpage()
      grid.draw(venn.plot)
      
      
      
    })
    
    
    
     # Render Virulence-Reistance heatmap plot based on the filter data sets
     output$Relative_abundance <- renderPlotly({
       tryCatch({
         # Merge abricate report and metadata 
         abri_kraken2_merged <- merge(filter_data(), metadata(), by ="SAMPLE")
         
       # Filter by Species and counts greater than 30
       df_species <- abri_kraken2_merged %>%
         group_by(NAME) %>%
         arrange(NAME) %>%
         group_by(SAMPLE, NAME, TREATMENT) %>%
         summarise(Count = n(), .groups = "drop") %>%
         group_by(SAMPLE) %>%
         mutate(Percentage = Count / sum(Count) * 100)
       
         # Annotation col colors for taxa
         # Description variable
         taxa <- as.factor(filter_data()$NAME)
         # print(taxa)
         
         # Generate a color palette based on the number of levels in gene_family
         t <- length(levels(taxa))
         print(t)
         t_palette <- paletteer_d("palettesForR::Cranes", n = t)
         t_palette
         
         
         
         #  stacked barplot where each bar represents a SAMPLE, 
         #  and the segments within each bar represent the relative abundance of different taxa
          ra <- ggplot(df_species, aes(x = factor(SAMPLE), y = Percentage, fill = NAME)) +
           geom_bar(stat = "identity") +
           scale_fill_manual(values = t_palette) +
           labs(x = "Sample groups", y = "Percentage", fill = "Taxonomic level") +
           theme_classic() +
           theme(legend.position="top") +
           theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
           theme(text = element_text(size = 14)) +
           facet_wrap(~ TREATMENT, scales = "free_x", nrow = 1)
         
         # Plot Coverage vs seq_length with PLOTLY ####
         Relative_abundace_plotly <- plotly::ggplotly(ra)
         Relative_abundace_plotly
         
         
         
       }, error = function(e) {
         # Handle the error
         # You can display an error message or take any other appropriate action
         showNotification("An error occurred while generating the plot. Number of requested colors greater than this palette can offer which is 256.",
                          type = "error")

         # Return an empty plot
         return(NULL)
       })
       
    })
     
     
     # Render ARG, VFs and MGEs distribution / Air and Surface SAMPLEs groups
     output$taxaCount_violinPlot <- renderPlot({
       tryCatch({
         # Merge abricate report and metadata 
         abri_kraken2_merged <- merge(dataframe(), metadata(), by ="SAMPLE")
         
         
         # Filter the data set by AMR, VFs and MGEs
         abri_kraken2_filtered <- abri_kraken2_merged %>%
           arrange(TREATMENT) %>%
           group_by(TREATMENT, NAME) %>%
           summarise(Taxa_count = n(), .groups = 'drop') %>%
           ungroup() %>% 
           mutate(Taxa_log_count = log(Taxa_count + 1))  # Adding 1 to avoid log(0)
         
         
         # one-way ANOVA analysis
         anova_result <- aov(Taxa_log_count ~ TREATMENT, data = abri_kraken2_filtered)
         summary(anova_result)
         
         # Generate the violin plot with statistical significance for taxa
         Taxa_count <- ggplot(abri_kraken2_filtered, aes(x = TREATMENT, y = Taxa_log_count, fill = TREATMENT)) +
           geom_violin(trim = FALSE, scale = "width") +
           geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
           scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                                        "Seaweed" = "#2ca02c",
                                        "Soyabean meal" = "#ff7f0e")) +
           labs(x = "Treatment groups", y = "Log(Number of Taxa)") +
           theme_classic() +
           theme(legend.position = "none", text = element_text(size = 16)) +
           annotate("text", x = 2, y = max(abri_kraken2_filtered$Taxa_log_count) + 1.3,
                    label = paste("ANOVA p =", format(summary(anova_result)[[1]][["Pr(>F)"]][1], digits = 2)),
                    size = 5, color = "black")+
           theme(text = element_text(size = 24))
         
         # Plot Taxa cont violin plot
         Taxa_count
         
         
       }, error = function(e) {
         # Handle the error
         # You can display an error message or take any other appropriate action
         showNotification("An error occurred while generating the plot. Adjust the level of each variable on the left to the minimum value.",
                          type = "error")

         # Return an empty plot
         return(NULL)
       })
       
       
     })
     
     # ARG Profile
     output$ARG_profile <- renderPlot({
       Gene_count <- table(filter_data()$DATABASE)
       
       # PLot the Gene count by database 
       pie_chart <- pie(Gene_count, labels = paste0(names(Gene_count), " / ", Gene_count, " / ", round(100 * Gene_count/sum(Gene_count), 2), "%"),
                        # main="Gene count by database",
                        col = c("#b02d02", "#5c0404", "#c47a02"),
                        cex = 1.3, # Adjust label size
                        radius = 1, # adjust the size
                        xlim = c(-1.5, 1.5) # Adjust x-axis limits to create more space for labels
       )
       
       # Plot pie chart
       pie_chart
                        
     })
     
     # ARG, VFs and MGEs by SAMPLEs
     output$GeneticsElements_SAMPLE <- renderPlot({
       # Count occurrences
       data_summary <- filter_data() %>%
         group_by(SAMPLE, DATABASE) %>%
         summarise(Count = n())
       
       # Calculate percentages
       data_summary <- data_summary %>%
         group_by(SAMPLE) %>%
         mutate(Percentage = Count / sum(Count) * 100)
       
       
       
       #  stacked barplot where each bar represents a SAMPLE, 
       #  and the segments within each bar represent the counts of different genes
       stacked_barplot <- ggplot(data_summary, aes(x = factor(SAMPLE), y = Percentage, fill = DATABASE)) +
         geom_bar(stat = "identity") +
         scale_fill_manual(values = c("#b02d02", "#5c0404", "#c47a02")) +
         labs(x = "Sample groups by Treatment", y = "Percentage", fill = "Database") +
         theme_classic() +
         theme(legend.position="top") +
         theme(text = element_text(size = 20)) +
         theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14)) +
         theme(axis.text.y = element_text(size = 14))
       
       # Plot the graph
       stacked_barplot
       
     })
    
     # Genetic elements co-occurrence on SAMPLEs groups
     output$geneticElements_occurrence <- renderPlot ({
       
       # Annotation rows
       # Pivot the SAMPLE values into columns
       df_wide_ann_rows <- dataframe() %>%
         arrange(SAMPLE) %>%
         pivot_wider(names_from = DATABASE, values_from = GENE, values_fn = list) %>%
         group_by(NAME) %>%
         summarise(across(15:17, ~ paste(na.omit(.), collapse = ", "))) %>%
         mutate(across(2:4, ~ str_replace_all(., "(NuLL,|,NULL|,NULL,|NULL|, )", "")))
       
       
       # Assign the values of NAMEs as row NAMEs
       df_wide_ann_rows <- df_wide_ann_rows %>%
         na.omit()  %>%
         remove_rownames %>%
         column_to_rownames(var="NAME") %>%
         mutate(across(everything(), ~ if_else(. == "", NA_character_, .)))
       
       # Convert from a dataframe to a matrix
       SAMPLE2arg_data_matrix <- data.matrix(df_wide_ann_rows, rownames.force = NA)
       
       # Count NA in data sets
       sum(is.na(SAMPLE2arg_data_matrix))
       
       # Replace all NA values with 0 values
       SAMPLE2arg_data_matrix <- SAMPLE2arg_data_matrix %>% replace(is.na(.), 0)
       
       # Convert to a binary matrix
       SAMPLE2arg_data_Bmatrix <- as.matrix((SAMPLE2arg_data_matrix > 0) + 0)
       
       # Create a list of sets
       sets <- apply(SAMPLE2arg_data_Bmatrix, 2, function(col) which(col == 1))
       names(sets) <- colnames(SAMPLE2arg_data_Bmatrix)
       print(sets)
       class(sets)
       
       ### Plot the Venn diagram
       venn.plot <- venn.diagram(
         x = sets,
         category.names = c("VFs", "AMR", "MGEs"),
         fill = c(vfdb="#EEAD0EFF", card="#CD3333FF", plasmidfinder="blue"), 
         alpha = 0.5,
         height = 50, width = 50,
         filename = NULL,
         cat.cex = 1.5,  # Adjust this value to change the font size
         cex = 2       # Adjust this value to change the font size of the numbers
       )
       grid.newpage()
       grid.draw(venn.plot)
       
       
     })
     
      
     # Render PCoA based on the ARG, VFs and MGEs profile
     output$PCoA_Analysis <- renderPlotly({
       suppressWarnings({
         
         # Merge abricate report and metadata 
         abri_kraken2_merged <- merge(dataframe(), metadata(), by ="SAMPLE")
         
         # Arrange dataset by gene variable
         abri_kraken2_arranged <- abri_kraken2_merged %>%
           arrange(GENE) %>%
           select(SAMPLE, NAME, GENE, TREATMENT)
         
         # Annotation rows
         # Pivot the Genes values into columns
         df_wide_ann_rows <- abri_kraken2_arranged %>%
           arrange(SAMPLE) %>%
           pivot_wider(names_from = SAMPLE, values_from = c(NAME), values_fn = length) %>%
           group_by(GENE) %>%
           summarise(across(2:19, ~ paste(., collapse = ", "))) %>%
           mutate(across(2:19, ~ str_replace_all(., "(NA,|,NA|,NA,|NA|,| )", "")))
         
         # Assign the values of GENE as row NAMEs and empty rows as NA value
         df_wide_ann_rows <- df_wide_ann_rows %>%
           na.omit()  %>%
           remove_rownames %>%
           column_to_rownames(var="GENE") %>%
           mutate(across(everything(), ~ if_else(. == "", NA_character_, .)))
         
         
         # Remove rows where all values are NA or empty strings
         df_clean <- df_wide_ann_rows[!apply(df_wide_ann_rows, 1, function(row) all(row == "" | is.na(row))), ]
         
         
         # Convert from a dataframe to a matrix
         SAMPLE2arg_data_matrix <- data.matrix(df_clean, rownames.force = NA)
         
         # Count NA in data sets
         sum(is.na(SAMPLE2arg_data_matrix))
         
         # Replace all NA values with 0 values
         SAMPLE2arg_data_matrix <- SAMPLE2arg_data_matrix %>% replace(is.na(.), 0)
         
         # Calculate Bray-Curtis distances
         bray_curtis_dist <- vegdist(SAMPLE2arg_data_matrix, method = "bray")
         
         # Perform PCoA
         pcoa_result <- cmdscale(bray_curtis_dist, eig = TRUE, k = 2)
         
         # Calculate percentage of variance explained
         eig_values <- pcoa_result$eig
         var_explained <- eig_values / sum(eig_values) * 100
         
         # Create a data frame for plotting
         pcoa_df <- data.frame(GENE = rownames(SAMPLE2arg_data_matrix),
                               PC1 = pcoa_result$points[, 1],
                               PC2 = pcoa_result$points[, 2])
         
         
         # Metadata dataframe
         TREATMENTData <- abri_kraken2_arranged %>%
           select(SAMPLE, GENE, TREATMENT) %>%
           distinct()
         
         # Merge the metadata and SAMPLE dataset
         # ann_pcoa <- inner_join(pcoa_df, MetadataLocations, by="SAMPLE")
         ann_pcoa <- inner_join(pcoa_df, TREATMENTData, by="GENE")
         
         # Annotation col colors for TREATMENT
         # Description variable
         TREATMENT_col <- as.factor(ann_pcoa$TREATMENT)
         
         # Generate a color palette based on the number of levels in gene_family
         l <- length(levels(TREATMENT_col))
         l_palette <- paletteer_d("ggsci::default_igv", n = l)
         l_palette
         
         # Plot the results
         pcaoa_plot <- ggplot(ann_pcoa, aes(x = PC1, y = PC2, colour = TREATMENT)) +
           geom_point(size = 3) +
           #geom_text(aes(label = GENE), vjust = 1.5, hjust = 1.5) +
           # stat_ellipse(geom = "polygon", aes(fill = TREATMENT), alpha = 0.2) +
           stat_ellipse(type = "norm", lwd = 0.8) + # Change line width
           # geom_jitter(width = 0.02, height = 0.02) +
           scale_color_manual(values = l_palette) +
           #scale_fill_manual(values = l_palette) +
           labs(title = "",
                color = "Treatment",
                x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
                y = paste0("PC2 (", round(var_explained[2], 2), "%)")
           ) +
           theme_classic() +
           theme(text = element_text(size = 16)) +
           scale_x_continuous(limits = c(-0.6, NA)) +  # Set x-axis to start at -0.5
           scale_y_continuous(limits = c(-0.8, NA))    # Set y-axis to start at -0.5
         
         # ggplot graph
         pcaoa_plot
         
         # ploty graph
         fig <- plotly::ggplotly(pcaoa_plot)
         fig
       
       })

     })
     
     # Render ARG, VFs and MGEs distribution / Air and Surface SAMPLEs groups
     output$ARGsCount_violinPlot <- renderPlot({
       tryCatch({
         # Merge abricate report and metadata 
         abri_kraken2_merged <- merge(dataframe(), metadata(), by ="SAMPLE")
         
         # Filter the data set by AMR
         abri_kraken2_filtered <- abri_kraken2_merged %>%
           #filter(grepl("card", DATABASE)) %>%
           arrange(NAME) %>%
           group_by(TREATMENT, GENE) %>%
           summarise(Gene_count = n(), .groups = 'drop') %>%
           ungroup() %>% 
           mutate(Gene_log_count = log(Gene_count + 1))  # Adding 1 to avoid log(0)
         
         
         # one-way ANOVA analysis
         anova_result <- aov(Gene_log_count ~ TREATMENT, data = abri_kraken2_filtered)
         summary(anova_result)
         
         
         # Generate the violin plot with statistical significance for AMR
         ggplot(abri_kraken2_filtered, aes(x = TREATMENT, y = Gene_log_count, fill = TREATMENT)) +
           geom_violin(trim = FALSE, scale = "width") +
           geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
           scale_fill_manual(values = c("Reference diet" = "#1f77b4",
                                        "Seaweed" = "#2ca02c",
                                        "Soyabean meal" = "#ff7f0e")) +
           labs(x = "Treatment groups", y = "Log(Number of GENE)") +
           theme_classic() +
           theme(legend.position = "none", text = element_text(size = 16),
                 axis.text.y = element_text(angle = 90, hjust = 0.5)
           ) +
           
           #geom_text(aes(label = Letters), vjust = -0.5, size = 5) +
           annotate("text", x = 2, y = max(abri_kraken2_filtered$Gene_log_count) + 0.3,
                    label = paste("ANOVA p =", format(summary(anova_result)[[1]][["Pr(>F)"]][1], digits = 2)),
                    size = 5, color = "black") +
           theme(text = element_text(size = 24))
         
         
         
       }, error = function(e) {
         # Handle the error
         # You can display an error message or take any other appropriate action
         showNotification("An error occurred while generating the plot. Adjust the level of each variable on the left to the minimum value.",
                          type = "error")
         
         # Return an empty plot
         return(NULL)
       })
       
       
     })
     
       #  # Key Pathogens-AMR co-occurrence
      output$Key_Pathoges_AMR <- renderPlot({
          # Merge bracken_arranged and metadata dataframe
        abri_kraken2_merged <- merge(dataframe(), metadata(), by ="SAMPLE")
        
        # Define modify_strings function
        modify_strings <- function(df, column_NAME, target_strings, replacement_strings) {
          if (length(target_strings) != length(replacement_strings)) {
            stop("target_strings and replacement_strings must be of the same length")
          }
          for (i in seq_along(target_strings)) {
            df <- df %>%
              mutate(!!sym(column_NAME) := str_replace_all(!!sym(column_NAME), target_strings[i], replacement_strings[i]))
          }
          return(df)
        }
        
        # List of target and replacement strings
        target_strings <- c("carbapenem;cephalosporin;penam")
        replacement_strings <- c("Beta-lactam")
        
        # Call Modify String function
        abri_kraken2_filtered <- modify_strings(abri_kraken2_merged, "RESISTANCE", target_strings, replacement_strings)
        
        # Filter and mutate the dataset
        abri_kraken2_filtered <- abri_kraken2_filtered %>%
          filter(grepl("card", DATABASE)) %>%
          mutate(GENE = case_when(
            grepl("vanW_gene_in_vanB_cluster", GENE) ~ "vanWB",
            grepl("vanR_gene_in_vanB_cluster", GENE) ~ "vanRB",
            grepl("vanX_gene_in_vanB_cluster", GENE) ~ "vanXB",
            grepl("vanY_gene_in_vanB_cluster", GENE) ~ "vanYB",
            grepl("vanS_gene_in_vanB_cluster", GENE) ~ "vanSB",
            grepl("vanH_gene_in_vanB_cluster", GENE) ~ "vanHB",
            grepl("Escherichia_coli_ampC_beta-lactamase", GENE) ~ "ampC",
            grepl("Escherichia_coli_mdfA", GENE) ~ "mdfA",
            grepl("Escherichia_coli_emrE", GENE) ~ "emrE",
            TRUE ~ GENE
          ),
          GENE = ifelse(grepl("AAC(6')-Ie-APH(2'')-Ia_bifunctional_protein", GENE), "aac(6')-Ie/aph(2'')-Ia", GENE)) %>%
          mutate(NAME = ifelse(grepl("Chryseobacterium panacisoli", NAME), "C. panacisoli", NAME)) %>%
          mutate(RESISTANCE = ifelse(nchar(RESISTANCE) > 20, "Multi-drug", RESISTANCE)) %>%
          arrange(NAME)
        
        # Clean the RESISTANCE strings
        abri_kraken2_filtered <- abri_kraken2_filtered %>%
          mutate(RESISTANCE = str_replace_all(RESISTANCE, "_", " ") %>% str_to_title())
        
        # Resistance gene category
        pathogen_AMR <- abri_kraken2_filtered %>%
          filter(!grepl("Alistipes|Azorhizobium|cellular organisms|Bacillota|Clostridia|Lachnospiraceae|Bacteria|root|Enterobacterales|Enterobacteriaceae|Terrabacteria group|Homo sapiens|Pseudomonadota|Bacteroidales|Bifidobacterium|Bacteroidota", NAME)) %>%
          filter(!(str_count(NAME, "\\S+") == 1 & NAME != "Enterococcus")) %>%
          filter(TREATMENT == "Reference diet") %>%
          # filter(TREATMENT == "Soyabean meal") %>%
          #filter(TREATMENT == "Seaweed") %>%
          select(SAMPLE, NAME, GENE, RESISTANCE) %>%
          mutate(NAME = str_replace(NAME, "^(\\w)\\w+\\s(\\w+)", "\\1. \\2")) %>%
          distinct()
        
        # Categories variables
        NAME_category <- unique(as.factor(pathogen_AMR$NAME))
        SAMPLE_category <- unique(as.factor(pathogen_AMR$SAMPLE))
        gene_category <- unique(as.factor(pathogen_AMR$GENE))
        resistant_category <- unique(as.factor(pathogen_AMR$RESISTANCE))
        
        #pathogen_AMR[26, 2] = "aac(6')-Ie/aph(2'')-Ia"
        # Reference diet
        pathogen_AMR[64, 2] = "C. innocuum"
        pathogen_AMR[65, 2] = "C. scindens"
        pathogen_AMR[67:68, 2] = "u. Subdoli."
        pathogen_AMR[66, 2] = "u. Plantact."
        
        # Filter by species
        # df <- abri_kraken2_filtered %>%
        #   filter(grepl("Clostridioides difficile|Enterococcus faecium|Staphylococcus aureus|Acinetobacter baumannii|Klebsiella pneumoniae|Klebsiella quasipneumoniae|Enterobacter cloacae complex|Enterobacter hormaechei|Escherichia coli|Streptococcus pneumoniae|Clostridium acetobutylicum|Clostridium perfringens|Neisseria cinerea|Streptococcus parasanguinis|Streptococcus pseudopneumoniae|Streptococcus sp. 116-D4|Streptococcus|Haemophilus|Moraxella|Staphylococcus|Corynebacterium|Neisseria|Prevotella|Veillonella", NAME))
        df <- pathogen_AMR
        
        # Create a connection matrix
        categories <- unique(c(sort(df$SAMPLE), sort(df$NAME), sort(df$GENE), sort(df$RESISTANCE)))
        mat <- matrix(0, nrow = length(categories), ncol = length(categories), dimnames = list(categories, categories))
        
        for (i in 1:nrow(df)) {
          mat[df$SAMPLE[i], df$NAME[i]] <- 1
          mat[df$NAME[i], df$GENE[i]] <- 1
          mat[df$GENE[i], df$RESISTANCE[i]] <- 1
        }
        
        # Set colors for categories
        color_sectors <- c(
          paletteer_d("ggsci::default_igv", n = length(unique(df$SAMPLE))),
          paletteer_d("ggsci::default_igv", n = length(unique(df$NAME))),
          paletteer_d("ggsci::default_igv", n = length(unique(df$GENE))),
          paletteer_d("ggsci::default_igv", n = length(unique(df$RESISTANCE)))
        )
        
        ## reset the graphic parameters and internal variables
        circos.clear()
        
        # Graph parameters
        circos.par(track.height = 0.1, 
                   start.degree = 115, 
                   gap.degree = 2, 
                   canvas.xlim = c(-1, 1), 
                   canvas.ylim = c(-1, 1), 
                   circle.margin = c(1, 1), 
                   unit.circle.segments = 500)
        
        # Create the chord diagram
        chordDiagram(mat, 
                     transparency = 0.5, 
                     annotationTrack = "grid", 
                     scale = FALSE, directional = 1, 
                     diffHeight = mm_h(3), 
                     grid.col = color_sectors, 
                     preAllocateTracks = list(track.height = 0.1, 
                                              unit.circle.segments = 500, 
                                              start.degree = 95, 
                                              scale = TRUE)
        )
        
        # Add labels to the sectors
        circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
          sector.index = get.cell.meta.data("sector.index")
          circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.8), cex = 1)
        }, bg.border = NA)
        
        
        
      })
      
      
      # Air microbiome heatmap
      output$Gut_microbiome <- renderPlot({
        # Merge bracken_arranged and metadata dataframe
        abri_kraken2_merged <- merge(dataframe(), metadata(), by ="SAMPLE")
        
        # Filter the dataset at the species levels
        abri_kraken2_clean <- abri_kraken2_merged %>%
          filter(!(str_count(NAME, "\\S+") == 1 & NAME != "Enterococcus"))
        
        # Pivot the SAMPLE values into columns
        df_wide <- abri_kraken2_clean %>%
          arrange(SAMPLE) %>%
          select(SAMPLE, TAXID, NAME, GENE, TREATMENT) %>%
          pivot_wider(names_from = SAMPLE, values_from = c(TAXID), values_fn = length) %>%
          group_by(NAME) %>%
          summarise(across(3:20, ~ paste(., collapse = ", "))) %>%
          mutate(across(2:19, ~ str_replace_all(., "(NA,|,NA|,NA,|NA|,| )", "")))
        
        # colnames(df_wide)
        
        # Assign the values of col1 as row NAMEs and handle empty strings
        processed_mgs2arg_data <- df_wide %>%
          remove_rownames() %>% 
          column_to_rownames(var = "NAME") %>%
          mutate(across(everything(), ~ na_if(., "")))
        
        # Remove rows where all values are NA
        SAMPLE2arg_data <- processed_mgs2arg_data %>%
          filter(rowSums(is.na(.)) < ncol(.))
        
        # Convert from a dataframe to a matrix and replace NA with 0
        SAMPLE2arg_data_matrix <- data.matrix(SAMPLE2arg_data, rownames.force = NA) %>%
          replace(is.na(.), 0)
        
        # Count NA in data sets
        sum(is.na(SAMPLE2arg_data_matrix))
        
        # Convert to a binary matrix
        SAMPLE2arg_data_Bmatrix <- as.matrix((SAMPLE2arg_data_matrix > 0) + 0)
        
        # Order the matrix by row NAMEs
        # SAMPLE2arg_data_Bmatrix <- SAMPLE2arg_data_Bmatrix[order(rownames(SAMPLE2arg_data_Bmatrix)), ]
        
        # Rearranging a Binary Matrix by Row Count
        rearrange_matrix <- function(mat) {
          # Calculate the row sums (count of '1's)
          ordered_indices <- order(rowSums(mat), decreasing = TRUE)
          
          # Rearrange the matrix
          return(mat[ordered_indices, ])
        }
        
        # Rearranging the binary matrix
        sorted_matrix <- rearrange_matrix(SAMPLE2arg_data_Bmatrix)
        
        # Modify ordering of the clusters using clustering callback option
        callback = function(hc, mat){
          sv = svd(t(mat))$v[,1]
          dend = reorder(as.dendrogram(hc), wts = sv)
          as.hclust(dend)
        }
        
        # Annotations col NAMEs
        # Transpose SAMPLE2arg_data
        SAMPLE2arg_data_transposed <- as.data.frame(t(SAMPLE2arg_data))
        
        # Convert rownames into column and select SAMPLE
        SAMPLE2arg_data_transposed <- SAMPLE2arg_data_transposed %>%
          mutate(SAMPLE = rownames(.)) %>%
          select(SAMPLE)
        
       # Select SAMPLE and treatmet category
        MetadataLocations <- abri_kraken2_clean %>%
          select(SAMPLE, TREATMENT) %>%
          distinct()
        
        # Merge the metadata and SAMPLE dataset
        ann_col <- inner_join(SAMPLE2arg_data_transposed, MetadataLocations, by = "SAMPLE") %>%
          na.omit() %>%
          remove_rownames() %>%
          column_to_rownames(var = "SAMPLE") 
        # rename(TYPE = type)
        
        # Generate a color palette based on the number of levels in location
        d_palette <- paletteer_d("ggsci::category10_d3", n = nlevels(as.factor(ann_col$TREATMENT)))
        
        # Create a NAMEd dataframe for the col colors
        df4 <- ann_col %>%
          distinct(TREATMENT) %>%
          mutate(color = d_palette)
        
        # Create a NAMEd vector for rows
        type_colors <- setNames(as.character(df4$color), df4$TREATMENT)
        
        # Filter the dataset by database
        # abri_kraken2_filtered <- abri_kraken2_merged[,c(-17,-18, -20, -21, -22)]
        
        # Annotation rows
        # Pivot the SAMPLE values into columns
        df_wide_ann_rows <- abri_kraken2_clean %>%
          arrange(SAMPLE) %>%
          filter(grepl("card|vfdb|plasmidfinder", DATABASE)) %>%
          pivot_wider(names_from = DATABASE, values_from = GENE, values_fn = list) %>%
          group_by(NAME) %>%
          summarise(across(16:18, ~ paste(., collapse = ", "))) %>%
          mutate(across(2:4, ~ str_replace_all(., "(NULL,|,NULL|,NULL,|NULL| )", ""))) %>%
          na.omit() %>%
          remove_rownames() %>%
          column_to_rownames(var = "NAME") %>%
          mutate(across(everything(), ~ if_else(. == "", NA_character_, .)))
        
        # Classify functions
        classify <- function(input_string, type) {
          if (is.na(input_string)) return(NA)
          elements <- unlist(strsplit(input_string, ","))
          unique_elements <- length(unique(elements))
          
          if (type == "virulence") {
            return(ifelse(unique_elements == 1, "Virulent", "Hypervirulent"))
          } else if (type == "amr") {
            return(ifelse(unique_elements == 1, "Drug-resistant", "Multidrug-resistant"))
          } else if (type == "MGE") {
            return(ifelse(unique_elements == 1, "Single MGE", "Multi-MGEs"))
          }
        }
        
        # Apply classification functions
        ann_rows <- df_wide_ann_rows %>%
          mutate(VFs = sapply(vfdb, classify, type = "virulence"),
                 AMR = sapply(card, classify, type = "amr"),
                 MGEs = sapply(plasmidfinder, classify, type = "MGE")) %>%
          select(VFs, AMR, MGEs)
        
        # Clean string function
        clean_string <- function(x) {
          str_to_title(str_replace_all(x, "_", " "))
        }
        
        # Apply the function to the RESISTANCE column
        Resistance_category <- abri_kraken2_clean %>%
          distinct(NAME, GENE, RESISTANCE) %>%
          mutate(RESISTANCE = sapply(RESISTANCE, clean_string)) %>%
          group_by(NAME) %>%  # Group by 'NAME' to join multiple values
          summarise(RESISTANCE = paste(unique(RESISTANCE), collapse = ", ")) %>%  # Join multiple values for RESISTANCE
          column_to_rownames(var = "NAME") %>%
          mutate(across(1, ~ str_replace_all(., "(NA,|,NA|,NA,|NA|,| )", ""))) %>%
          mutate(RESISTANCE = ifelse(nchar(RESISTANCE) > 20, "Multi-drug", RESISTANCE)) %>%
          mutate(across(everything(), ~ if_else(. == "", NA_character_, .))) %>%
          select(RESISTANCE)
        
        # Merge all the annotation rows
        ann_row_all <- merge(ann_rows, Resistance_category, by = 0) %>%
          rename(DRUG = RESISTANCE) %>%
          column_to_rownames(var = "Row.names") %>%
          select(DRUG, AMR, VFs, MGEs)
        
        # Generate a color palette for drug categories
        r_palette <- paletteer_d("ggsci::default_igv", n = length(unique(ann_row_all$DRUG))-1)
        
        # Create a NAMEd dataframe for the Resistance category
        ann_row_drug_col <- ann_row_all %>%
          select(DRUG) %>%
          na.omit() %>%
          distinct(DRUG) %>%
          mutate(color = r_palette)
        
        # Create a NAMEd vector for rows
        drug_colors <- setNames(as.character(ann_row_drug_col$color), ann_row_drug_col$DRUG)
        
        # Annotation colors customization
        ann_colors <- list(
          TREATMENT = type_colors,
          DRUG = drug_colors,
          AMR = c(`Drug-resistant` = "#a02d02", `Multidrug-resistant` = "#522501"),
          MGEs = c(`Multi-MGEs` = "#5c0404", `Single MGE` = "#fcd2d2"),
          VFs = c(`Virulent` = "#fce5d2", `Hypervirulent` = "#c47a02")
        )
        
        
        # Create heatmap using pheatmap package ##
        # heatmap_plot <- pheatmap(sorted_matrix[1:50,], display_numbers = FALSE, cluster_cols = TRUE, cluster_rows = FALSE,
        #                          scale = "none", 
        #                          clustering_callback = callback,  
        #                          border_color = "NA", color = c("#CCCCCCFF", "#666666FF"),
        #                          legend_breaks = c(0, 1),
        #                          legend_labels = c("Absent", "Present"),
        #                          annotation_row = ann_row_all,
        #                          annotation_col = ann_col[1],
        #                          show_rownames = TRUE,
        #                          # cutree_cols = 15,
        #                          # cutree_rows = 20,
        #                          annotation_colors = ann_colors,
        #                          fontsize_row = 14,  # Adjust this value for row NAMEs
        #                          fontsize_col = 14,   # Adjust this value for column NAMEs
        #                          fontsize = 14  # Adjust this value for the legend text
        #                          
        #                          
        # )
        # 
        # # Plot gut microbiome pheatmap
        # heatmap_plot
        
       # Attempt to create the heatmap
       tryCatch({
         pheatmap(sorted_matrix[input$row_count[1]:input$row_count[2], ], display_numbers = FALSE, cluster_cols = TRUE, cluster_rows = FALSE,
                  scale = "none", 
                  clustering_callback = callback,
                  border_color = "NA", color = c("#CCCCCCFF", "#666666FF"),
                  legend_breaks = c(0, 1),
                  legend_labels = c("Absent", "Present"),
                  annotation_row = ann_row_all,
                  annotation_col = ann_col[1],
                  show_rownames = TRUE,
                  # cutree_cols = 15,
                  # cutree_rows = 20,
                  annotation_colors = ann_colors,
                  fontsize_row = 14,  # Adjust this value for row NAMEs
                  fontsize_col = 14,   # Adjust this value for column NAMEs
                  fontsize = 14  # Adjust this value for the legend text


         )
       }, error = function(e) {
         # Handle error: display a message in the console
         message("An error occurred while generating the heatmap: ", e$message)

       })
       
       
        
        
      })

    
    # Return the reactive that yields the data frame
    return(filter_data)
  }
  
  )
}
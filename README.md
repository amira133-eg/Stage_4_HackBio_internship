# Gene Ontology Enrichment Analysis Shiny App
م
## Overview
This Shiny app performs **Gene Ontology (GO) enrichment analysis** for a set of genes known to be upregulated in cancer. It allows users to either use a predefined set of genes or upload their own gene list. The app utilizes the **TCGAanalyze_EAcomplete()** and **TCGAvisualize_EAbarplot()** functions from the `TCGAbiolinks` package to perform the analysis and display a bar plot of enriched GO terms


app. link:[ Enrichment Analysis](https://enrichmentanalysis.shinyapps.io/TCGAapp/)
---

## Instructions on How to Use the App

1. **Predefined Gene List**:
   - The app provides a predefined list of genes that are commonly associated with cancer.

2. **Upload My Genes**:
   - Users can input their own gene list, either in Gene Symbol or Ensembl ID format.

3. **Ontology Categories**:
   - Users can choose from three GO categories: Biological Process (BP), Cellular Component (CC), and Molecular Function (MF).

4. **Output**:
   - The output includes a bar plot displaying the top 10 enriched GO terms for the selected ontology.

---

				
# Features

- **Predefined Gene List**: Users can analyze a list of genes known to be upregulated in cancer.
- **Custom Gene List**: Users can upload their own gene list (in either gene symbol or Ensembl ID format).
- **Ontology Categories**: Choose from Biological Process, Cellular Component, or Molecular Function.
- **Visualization**: Results are displayed as a bar plot showing the top 10 enriched GO terms for the selected ontology.


			
## Code Explanation
			
# Load necessary libraries for building the Shiny app and performing gene analysis
      library(shiny)        # Main library for creating Shiny web applications
      library(TCGAbiolinks) # Library to perform cancer genomics analysis using TCGA data
      library(org.Hs.eg.db) # Annotation package for human gene data
      library(AnnotationDbi)# Interface for querying gene annotations
     library(rsconnect)    # Library for deploying Shiny applications
     library(BiocManager)  # Bioconductor package management

# List of predefined genes (gene symbols) commonly associated with cancer
     predefined_genes <- c(
      "TP53", "BRCA1", "BRCA2", "PIK3CA", "PTEN", "CDH1", "AKT1", "EGFR", "RB1", "ERBB2",
      "CCND1", "CDK4", "VEGFA", "ESR1", "MYC", "MDM2", "AR", "KRAS", "NRAS", "JUN",
    # (additional genes in the list) ...
     "HSP90AA1", "HSPB1", "GSTM1", "TP73"
     )

# Function to convert gene symbols (or other types of IDs) into Entrez IDs using org.Hs.eg.db
    convert_to_entrez <- function(gene_list, id_type = "SYMBOL") {
     entrez_ids <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = gene_list, column = "ENTREZID", 
    keytype = id_type, multiVals = "first"
     )
        entrez_ids <- na.omit(entrez_ids)  # Remove any entries with NA values
     as.character(entrez_ids)           # Return Entrez IDs as character vectors
    }

# Define the UI (user interface) for the Shiny app
     ui <- fluidPage(
    titlePanel("Gene Ontology Enrichment Analysis"),  # Title of the app
  
     # Tabs for navigating different parts of the app
    tabsetPanel(
    # Instructions tab explaining how to use the app
    tabPanel("Instructions",
             h3("How to Use This App"),
             p("This app performs Gene Ontology (GO) enrichment analysis for a set of genes..."),
             h4("Steps to Use"),
             p("1. In the 'Analysis' tab, select..."),
             tags$ul(
               tags$li("Enter your gene list in a comma-separated format..."),
               # More details on how to enter and format genes
             )
    ),
    
    # Analysis tab where users input gene lists and run analysis
    tabPanel("Analysis",
             sidebarLayout(
               sidebarPanel(
                 h4("Choose Analysis Parameters"),
                 # Radio buttons to choose between predefined or custom gene list
                 radioButtons("geneSource", "Gene List Source", choices = c("Predefined Gene List", "Upload My Genes")),
                 # Conditional input for custom gene list and type of gene IDs
                 conditionalPanel(
                   condition = "input.geneSource == 'Upload My Genes'",
                   textAreaInput("customGenes", "Enter Your Gene List (comma-separated)", placeholder = "e.g., BRCA1, TP53, MYC"),
                   selectInput("idType", "ID Type", choices = c("Gene Symbol" = "SYMBOL", "Ensembl ID" = "ENSEMBL"))
                 ),
                 # Dropdown to choose ontology type
                 selectInput("ontology", "Choose Ontology", choices = c("Biological Process" = "BP", "Cellular Component" = "CC", "Molecular Function" = "MF")),
                 actionButton("analyze", "Run Enrichment Analysis"),  # Button to run analysis
                 br()
               ),
               mainPanel(
                 h4("Note"),
                 p("After clicking 'Run Enrichment Analysis', the bar plot will display in the 'Results' tab.")
               )
             )
    ),
    
    # Results tab to display the output of the enrichment analysis as a bar plot
    tabPanel("Results",
             h4("Enrichment Analysis Results"),
             plotOutput("barPlot")  # Area for showing the bar plot
    )
    )
    ) 
 
# Server logic to handle user inputs and perform the analysis
    server <- function(input, output) {
  
    # Event listener for the 'analyze' button
    observeEvent(input$analyze, {
    # Decide whether to use predefined or custom genes
    gene_list <- if (input$geneSource == "Predefined Gene List") {
      predefined_genes
    } else {
      strsplit(input$customGenes, ",\\s*")[[1]]  # Split custom gene list by commas
    }
    
    # Convert the gene list to Entrez IDs
    entrez_ids <- convert_to_entrez(gene_list, id_type = input$idType)
    
    # If no valid Entrez IDs, show a message
    if (length(entrez_ids) == 0) {
      output$barPlot <- renderPlot({
        plot.new()
        title("No valid Entrez IDs found. Please check your input.")
      })
      return()
    }
    
    # Perform enrichment analysis using the TCGA functions
    enrichment_results <- TCGAanalyze_EAcomplete(
      TFname = "Breast Cancer Analysis",
      RegulonList = list(entrez_ids)
    )
    
    # Display the enrichment results in a bar plot
    output$barPlot <- renderPlot({
      ontology_result <- switch(input$ontology,
                                BP = enrichment_results$ResBP,
                                CC = enrichment_results$ResCC,
                                MF = enrichment_results$ResMF)
      
      if (!is.null(ontology_result) && nrow(ontology_result) > 0 && !any(is.na(ontology_result))) {
        # Visualize the top 10 GO terms
        TCGAvisualize_EAbarplot(GOBPTab = ontology_result, nBar = 10, color = "skyblue")
      } else {
        plot.new()
        title("No Enrichment Results Found")
      }
    })
     })
     }

# Run the app
    shinyApp(ui = ui, server = server)


               

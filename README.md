# Gene Ontology Enrichment Analysis Shiny App
# **Authors** 
(@slack):Amira Mahmoud Mohamed (amira_mahmoud_4463)

Manar Tarek (ManarTj)

Noran Morad (NoranMorad)


app. link:[ Enrichment Analysis](https://enrichmentanalysis.shinyapps.io/enrichmentgo/)
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
     library(shiny)           # For creating the interactive web application
     library(org.Hs.eg.db)     # Human gene annotations, including Entrez IDs and gene symbols
     library(AnnotationDbi)    # Provides functions for mapping gene IDs (e.g., SYMBOL to ENTREZID)
     library(clusterProfiler)  # For gene set enrichment analysis, including GO and KEGG pathways
     library(enrichplot)       # For visualizing the results of enrichment analysis
     library(DT)              # For creating interactive tables (DataTables)


# List of predefined genes (gene symbols) commonly associated with cancer
       example_gene_list <- c("TP53", "BRCA1", "BRCA2", "PIK3CA", "PTEN", "CDH1", "AKT1", 
                       "EGFR", "RB1", "ERBB2", "CCND1", "CDK4", "VEGFA", "ESR1", 
                       "MYC", "MDM2", "AR", "KRAS", "NRAS", "JUN")

# Function to convert gene symbols to Entrez IDs (used in enrichment analysis)
# gene_list: List of gene symbols or IDs
# id_type: Type of ID provided (e.g., SYMBOL or ENSEMBL)
     convert_to_entrez <- function(gene_list, id_type = "SYMBOL") {
     if (length(gene_list) > 0) {
    entrez_ids <- AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = gene_list, column = "ENTREZID", 
      keytype = id_type, multiVals = "first"  # Map each gene symbol to its Entrez ID
    )
    na.omit(as.character(entrez_ids))  # Remove any NA entries (unmapped genes)
    } else {
    character(0)  # Return empty if no genes are provided
    }
    }


# Define the UI (user interface) for the Shiny app
       ui <- fluidPage(
     titlePanel("Gene Ontology and Pathway Enrichment Analysis"),  # App title
  
     # Tabs to organize content: Guide, Analysis, and Results
     tabsetPanel(
    tabPanel("Guide",  # First tab explains how to use the app
             h3("How to Use This App"),
             p("This app performs Gene Ontology (GO) and KEGG Pathway enrichment analysis for a set of genes."),
             p("Users can analyze the predefined gene list or enter their own list."),
             h4("Steps to Use"),  # Step-by-step instructions for users
             tags$ol(
               tags$li("Go to the 'Analysis' tab and select 'Example Gene List' or 'Custom Gene List'."),
               tags$li("If using a custom list, enter your gene symbols or Ensembl IDs, separated by commas."),
               tags$li("Select the type of ID (Gene Symbol or Ensembl ID) from the dropdown."),
               tags$li("Click 'Run Enrichment Analysis' to see the results in the 'Results' tab.")
             ),
             h4("Expected Output"),  # Description of the expected results
             p("The results include bar plots for the top 10 enriched GO terms in each ontology: Biological Process, Molecular Function, and Cellular Component."),
             p("Additionally, KEGG pathway enrichment results are visualized and provided in interactive tables.")
    ),
    
    # Tab for performing the analysis
    tabPanel("Analysis",
             sidebarLayout(
               sidebarPanel(
                 h4("Parameters"),
                 radioButtons("geneSource", "Gene List Source", choices = c("Example Gene List", "Custom Gene List")),  # Choose between example or custom genes
                 conditionalPanel(
                   condition = "input.geneSource == 'Custom Gene List'",
                   textAreaInput("customGenes", "Enter Your Gene List (comma-separated)", 
                                 placeholder = "e.g., TP53, BRCA1, MYC"),  # Input box for custom gene list
                   selectInput("idType", "ID Type", choices = c("Gene Symbol" = "SYMBOL", "Ensembl ID" = "ENSEMBL"))  # Select the ID type
                 ),
                 actionButton("analyze", "Run Enrichment Analysis"),  # Button to start the analysis
                 br(), br(),
                 p("Note: Enrichment analysis may take a few seconds to complete.")
               ),
               mainPanel(
                 h4("Note"),
                 p("The results will be displayed in the 'Results' tab once analysis is complete.")  # Instructional text
               )
             )
    ),
    
    # Results tab with multiple sub-tabs for visualizing GO and KEGG results
    tabPanel("Results",
             h4("Enrichment Analysis Results"),
             tabsetPanel(
               tabPanel("GO Plots",  # Sub-tab for GO enrichment plots
                        h5("Biological Process"),
                        plotOutput("plot_BP"),  # Plot for biological process
                        h5("Molecular Function"),
                        plotOutput("plot_MF"),  # Plot for molecular function
                        h5("Cellular Component"),
                        plotOutput("plot_CC")   # Plot for cellular component
               ),
               tabPanel("Pathway Plots",  # Sub-tab for KEGG pathway plots
                        h5("KEGG Pathway Enrichment"),
                        plotOutput("plot_KEGG")  # Plot for KEGG pathway enrichment
               ),
               tabPanel("GO Tables",  # Sub-tab for GO enrichment tables
                        h5("Biological Process"),
                        DTOutput("table_BP"),   # Table for biological process results
                        h5("Molecular Function"),
                        DTOutput("table_MF"),   # Table for molecular function results
                        h5("Cellular Component"),
                        DTOutput("table_CC")    # Table for cellular component results
               ),
               tabPanel("Pathway Table",  # Sub-tab for KEGG pathway tables
                        h5("KEGG Pathway Enrichment"),
                        DTOutput("table_KEGG")  # Table for KEGG pathway results
               )
             )
    )
    )
    )

 
# Server logic to handle user inputs and perform the analysis
        server <- function(input, output) {
     # Create reactive values to store enrichment results
    enrichment_results <- reactiveValues(
    enriched_go_BP = NULL,
    enriched_go_MF = NULL,
    enriched_go_CC = NULL,
    enriched_KEGG = NULL
    )
  
    # Observe the event when the "analyze" button is clicked
    observeEvent(input$analyze, {
    
    # Get the gene list based on user input
    gene_list <- if (input$geneSource == "Example Gene List") {
      example_gene_list  # Use the predefined gene list
    } else {
      # Get the user-provided gene list, split it by commas, and clean it up
      genes_input <- unlist(strsplit(input$customGenes, ",\\s*"))
      genes_input[genes_input != ""]  # Remove any empty entries
    }
    
    # Convert the gene list to Entrez IDs (needed for enrichment analysis)
    entrez_ids <- convert_to_entrez(gene_list, id_type = input$idType)
    
    # Show a message if no valid Entrez IDs are found
    if (length(entrez_ids) == 0) {
      showModal(modalDialog(
        title = "No valid Entrez IDs found",
        "Please check your input and ensure you selected the correct ID type."
      ))
      return()
    }
    
    # Perform GO enrichment analysis for each ontology (BP, MF, CC)
    enrichment_results$enriched_go_BP <- enrichGO(gene = entrez_ids, 
                                                  OrgDb = org.Hs.eg.db, 
                                                  ont = "BP",  # Biological Process
                                                  pAdjustMethod = "BH",  # Adjust p-values
                                                  pvalueCutoff = 0.05, 
                                                  qvalueCutoff = 0.2, 
                                                  readable = TRUE)  # Convert IDs to gene symbols for readability
    enrichment_results$enriched_go_MF <- enrichGO(gene = entrez_ids, 
                                                  OrgDb = org.Hs.eg.db, 
                                                  ont = "MF",  # Molecular Function
                                                  pAdjustMethod = "BH", 
                                                  pvalueCutoff = 0.05, 
                                                  qvalueCutoff = 0.2, 
                                                  readable = TRUE)
    enrichment_results$enriched_go_CC <- enrichGO(gene = entrez_ids, 
                                                  OrgDb = org.Hs.eg.db, 
                                                  ont = "CC",  # Cellular Component
                                                  pAdjustMethod = "BH", 
                                                  pvalueCutoff = 0.05, 
                                                  qvalueCutoff = 0.2, 
                                                  readable = TRUE)
    
    # Perform KEGG pathway enrichment analysis
    enrichment_results$enriched_KEGG <- enrichKEGG(gene = entrez_ids, 
                                                   organism = 'hsa',  # Human organism ('hsa')
                                                   pAdjustMethod = "BH", 
                                                   pvalueCutoff = 0.05, 
                                                   qvalueCutoff = 0.2)
    
     })  # Close the observeEvent for the analyze button
  
    # Render the GO Biological Process plot
     output$plot_BP <- renderPlot({
     req(enrichment_results$enriched_go_BP)  # Wait for the results
     if (nrow(enrichment_results$enriched_go_BP@result) > 0) {
       barplot(enrichment_results$enriched_go_BP, showCategory = 10, title = "Biological Process")  # Plot the top 10 GO terms
    }
     })
  
    # Render the GO Molecular Function plot
     output$plot_MF <- renderPlot({
    req(enrichment_results$enriched_go_MF)
    if (nrow(enrichment_results$enriched_go_MF@result) > 0) {
      barplot(enrichment_results$enriched_go_MF, showCategory = 10, title = "Molecular Function")
    }
    })
  
    # Render the GO Cellular Component plot
     output$plot_CC <- renderPlot({
    req(enrichment_results$enriched_go_CC)
    if (nrow(enrichment_results$enriched_go_CC@result) > 0) {
      barplot(enrichment_results$enriched_go_CC, showCategory = 10, title = "Cellular Component")
    }
    })
  
    # Render the KEGG Pathway plot
    output$plot_KEGG <- renderPlot({
    req(enrichment_results$enriched_KEGG)
    if (nrow(enrichment_results$enriched_KEGG@result) > 0) {
      barplot(enrichment_results$enriched_KEGG, showCategory = 10, title = "KEGG Pathway Enrichment")
    }
    })
  
     # Render the GO Biological Process table
     output$table_BP <- renderDT({
    req(enrichment_results$enriched_go_BP)
    datatable(as.data.frame(enrichment_results$enriched_go_BP@result), options = list(pageLength = 5, autoWidth = TRUE))
     })
  
    # Render the GO Molecular Function table
     req(enrichment_results$enriched_go_MF)
    datatable(as.data.frame(enrichment_results$enriched_go_MF@result), options = list(pageLength = 5, autoWidth = TRUE))
    })
   
    # Render the GO Cellular Component table
    output$table_CC <- renderDT({
    req(enrichment_results$enriched_go_CC)
    datatable(as.data.frame(enrichment_results$enriched_go_CC@result), options = list(pageLength = 5, autoWidth = TRUE))
    })
  
    # Render the KEGG Pathway table
    output$table_KEGG <- renderDT({
    req(enrichment_results$enriched_KEGG)
    datatable(as.data.frame(enrichment_results$enriched_KEGG@result), options = list(pageLength = 5, autoWidth = TRUE))
    })
  
    }

# Run the app
    shinyApp(ui = ui, server = server)


               

# martGO2list.shinyapp# A R/shiny tool to query BioMart with a GO term and get a gene list# the ENSGID list can be used to plot a heatmap with fpkm2heatmap# SP@NC, 2024-03-08

library(shiny)
library(shinyBS)
library(shinyjs)
library(biomaRt)
library(DT)
library(openxlsx)

# user inputs: organism and GO:term 
# search for all protein_coding genes linked to the GO term
# output a table with information
# export ENSGID list to text file for use in fpkm1heatmap
# export full table to excel if needed

# you may uncomment the next line to allow large input files
options(shiny.maxRequestSize=1000*1024^2)
# the following test checks if we are running on shinnyapps.io to limit file size dynamically
# ref: https://stackoverflow.com/questions/31423144/how-to-know-if-the-app-is-running-at-local-or-on-server-r-shiny/31425801#31425801
#if ( Sys.getenv('SHINY_PORT') == "" ) { options(shiny.maxRequestSize=1000*1024^2) }

app.name <- "martGO2list"
script.version <- "1.0.0"

# Define UI
ui <- fluidPage(
  useShinyjs(),
  HTML('<style type="text/css">
    .row-fluid { width: 20%; }
    .well { background-color: #99CCFF; }
    .shiny-html-output { font-size: 14px; line-height: 15px; }
    .disabled-button {
       pointer-events: none; /* Disables click events */
       color: #ccc; /* Grey out the button text */
       background-color: #eee; /* Grey out the button background */
       }
       </style>'),
  # Application header
  headerPanel("Query the BioMart gene DB with a GO term"),

  # Application title
  titlePanel(
    windowTitle = "BioMart GO Data Retrieval",
    tags$a(href="https://nucleomicscore.sites.vib.be/en", target="_blank",
           img(src='logo.png', align = "right", width="150", height="58.5", alt="VIB Nucleomics Core"))
  ),

  # Sidebar with input
  sidebarLayout(
    # show file import and molecule filters
    sidebarPanel(
      tags$h5(paste(app.name, " version: ", script.version, sep="")),
      tipify(selectInput("enshost", "Select a host:",
              choices = c("https://www.ensembl.org",
                          "https://useast.ensembl.org",
                          "https://asia.ensembl.org"),
              selected="https://useast.ensembl.org"), # Default to EU main site
          "if the Main site appears Offline, try a mirror"),
      tipify(selectInput("organism", "Select an organism:",
                          choices = NULL), # This will be updated dynamically
             "all current species are shown in the list"),
      tipify(selectInput("biotypes", "Select one or more biotypes (click to add, Del key to remove):",
                  choices = NULL, multiple = TRUE), # This will be updated dynamically
             "by default only coding genes are selected"),
      tipify(textInput("go_term", "Enter GO Parent Term:", placeholder = "e.g., GO:0008150"),
             "use QuickGO to identify a relevant GO term"),
      actionButton("run_query", "Run"),
      tags$br(),
      tags$br(),
      downloadButton("save_list_button", "Save Ensembl IDs", class = "disabled-button"),
      downloadButton("save_table_button", "Save Table (txt)", class = "disabled-button"),
      downloadButton("save_table_button2", "Save Table (xlsx)", class = "disabled-button"),
      tags$br(),
      tags$br(),
      tipify(downloadButton('downloadMM', 'Download M&M'),
             "Download a text including the names and versions of all packages used in this webtool")
    ),
    # Show a plot of the generated distribution
    mainPanel(
      textOutput("selected_biotypes"),
      textOutput("row_count"),
      tags$br(),
      DTOutput("data_table")
    )
  )
# end UI bloc
)

# Define server
server <- function(input, output, session) {
  
  # Observe the data table and enable/disable download buttons accordingly
  observe({
    if (nrow(gene_data()) > 0) {
      # Enable the buttons if the table has data
      shinyjs::removeClass(selector = "#save_list_button", class = "disabled-button")
      shinyjs::removeClass(selector = "#save_table_button", class = "disabled-button")
      shinyjs::removeClass(selector = "#save_table_button2", class = "disabled-button")
    } else {
      # Disable the buttons if the table is empty
      shinyjs::addClass(selector = "#save_list_button", class = "disabled-button")
      shinyjs::addClass(selector = "#save_table_button", class = "disabled-button")
      shinyjs::addClass(selector = "#save_table_button2", class = "disabled-button")
    }
  })
  
  # Observe the selected host and update the organism choices
  observe({
    updateSelectInput(session, "organism",
                      choices = listDatasets(useMart(biomart = "ensembl", host=input$enshost))$dataset,
                      selected = "hsapiens_gene_ensembl")  # Default to Human)
  })
  
  # Reactive expression to update biotypes based on the selected organism
  observeEvent(input$organism, {
    if (!is.null(input$organism)) {
      # Ensure the dataset is valid before calling keys
      valid_datasets <- listDatasets(useMart(biomart = "ensembl", host=input$enshost))$dataset
      if (input$organism %in% valid_datasets) {
        updateSelectInput(session, "biotypes",
                          choices = as.list(keys(useMart(biomart = "ensembl", host=input$enshost, dataset=input$organism), keytype="biotype")),
                          selected = "protein_coding")  # Default to protein_coding
      }
    }
  }, ignoreNULL = TRUE)
  
  # create connection to Ensembl
  ensembl_mart <- eventReactive(input$run_query, {
    invalidateLater(1000)  # Simulate delay for demonstration
    useMart(biomart = "ensembl", host=input$enshost, dataset = input$organism)
  })
  
  # list biotypes
  output$selected_biotypes <- renderText({
    paste("Selected biotypes:", paste(input$biotypes, collapse = ", "))
  })
  
  # Retrieve data based on GO term and biotype "protein_coding"
  gene_data <- eventReactive(input$run_query, {
    if (!is.null(input$go_term)) {
      query <- input$go_term
      biotype <- as.vector(input$biotypes)
      getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype", 
                           "chromosome_name", "start_position",
                           "end_position", "strand"),
            filters = c("go_parent_term", "biotype"),
            values = list(query, biotype),
            mart = ensembl_mart())
    } else {
      data.frame()  # Empty data frame if no GO term provided
    }
  })
  
  # Display data in a DataTable (excluding biotype)
  output$data_table <- renderDT({
    datatable(gene_data(), rownames = FALSE)
  })
  
  # Render the number of rows in the data table
  output$row_count <- renderText({
    paste("Number of rows in table:", nrow(gene_data()))
  })
  
  # Save Ensembl IDs to a text file
  output$save_list_button <- downloadHandler(
    filename = function() {
      paste0(gsub(":", "_", input$go_term), "_", sub("_(.*)", "", input$organism), "_list.txt")
    },
    content = function(file) {
      writeLines(gene_data()$ensembl_gene_id, con = file)
    }
  )
  
  # Save table to an Excel file
  output$save_table_button <- downloadHandler(
    filename = function() { 
      paste0(gsub(":", "_", input$go_term), "_", sub("_(.*)", "", input$organism), "_info.txt", sep="")
    },
    content = function(file) {
      write.table(as.data.frame(gene_data()), file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
    }
  )
  
  # Save table to an Excel file
  output$save_table_button2 <- downloadHandler(
    filename = function() { 
      paste0(gsub(":", "_", input$go_term), "_", sub("_(.*)", "", input$organism), "_info.xlsx", sep="")
    },
    content = function(file) {
      write.xlsx(as.data.frame(gene_data()), file, colNames = TRUE, asTable = TRUE, overwrite = TRUE, colWidths = "auto", withFilter = TRUE)
    }
  ) 

  # download M&M text
  output$downloadMM <- downloadHandler(
    filename = function() {
      paste(app.name, "session_info.txt", sep="_")
    },
    content = function(file) {
      sink(file, append=TRUE)
      cat(paste("Thanks for using our tool", app.name, script.version, "\n", sep=" "))
      cat ("This tool sends a query to BioMart for a GO:term and retrieve information about the corresponding genes.")
      cat("\nYou can contact the Nucleomics Core at nucleomics.bioinformatics@vib.be with your question\n")
      cat(paste("This data was generated on ", format(Sys.time(), "%a %b %d %H:%M:%S %Y"), "\n",sep=" "))
      cat("\nthe R packages used in the tools are listed next:\n")
      print(capture.output(sessionInfo()))
      sink()
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)

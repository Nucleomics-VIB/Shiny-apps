# RNASeqFiltering.shinyapp
# A R/shiny tool to filter RNASeq statistical analysis results
# This tool outputs a list of genes (signature)
# that can be used with the companion tool fpkm2heatmap
# together with the matching Nucleomics Core MS-Excel RNASeq count file

library("shiny")
library("shinyBS")
library("openxlsx")
library("DT")
library("pheatmap")
library("RColorBrewer")

# you may uncomment the next line to allow large input files
options(shiny.maxRequestSize=1000*1024^2)
# the following test checks if we are running on shinnyapps.io to limit file size dynamically
# ref: https://stackoverflow.com/questions/31423144/how-to-know-if-the-app-is-running-at-local-or-on-server-r-shiny/31425801#31425801
#if ( Sys.getenv('SHINY_PORT') == "" ) { options(shiny.maxRequestSize=1000*1024^2) }

app.name <- "RNASeqFiltering"
script.version <- "1.0"

#defaults
def.min.lfc <- 0
def.max.pv <- 1
contrast.names.all <- NULL

# Define UI for application that draws a histogram
ui <- fluidPage(
  HTML('<style type="text/css">
       .row-fluid { width: 25%; }
       .well { background-color: #99CCFF; }
       .shiny-html-output { font-size: 14px; line-height: 15px; }
       </style>'),
  # Application header
  headerPanel(tags$h3("Filter RNASeq statistical data based on LogFC and corPvalue limits")),

  # Application title
  titlePanel(
    windowTitle = "RNASeq Filtering Tool",
    tags$a(href="https://corefacilities.vib.be/nc", target="_blank",
           img(src='logo.png', align = "right", width="150", height="58.5", alt="VIB Nucleomics Core"))
  ),

  # Sidebar with input
  sidebarLayout(
    # show file import and molecule filters
    sidebarPanel(
      tags$h5(paste(app.name, " version: ", script.version, sep="")),
      downloadButton("downloadData", label = "Download test data"),
      tags$br(),
      tags$a(href="license.pdf", target="_blank", "usage licence"),
      tags$hr(),
      tipify(fileInput('file1', 'Choose a StatisticalResults XLSX File', accept='.xlsx'),
             "the Data is a MS-Excel file provided by the Nucleomics Core, with a single worksheet reporting pairwise group comparison results (with logFC and FDR columns for each contrast), you may produce a compatible file based on the test data provided here."),
      tags$h4("modify cutoffs and click ", tags$em("Filter")),
      tipify(uiOutput("choose_contrasts"),"Use the checkboxes to define which contrasts should be filtered to produce the result table (all by default)"),
      tipify(textInput('min.lfc', "logFC cutoff", value = def.min.lfc, width = "150px"),"provide the absolute value of the logFC cutoff. All genes whose logFC is < -(value) or > +(value) will be kept"),
      tipify(textInput('max.pv', "corPvalue limit", value = def.max.pv, width = "150px"),"provide the max acepted FDR (corrected Pvalue) to keep a gene"),
      tipify(radioButtons('stringency', "Filter stringency:", c("one_or_more", "all_contrasts"), selected="one_or_more"),"the filtering should apply to at least one of the contrasts or to all selected contrasts (default to at least one)"),
      actionButton(inputId='goButton', "Filter", style='padding:4px; font-weight: bold; font-size:150%'),
      tags$hr(),
      tags$a(href="javascript:history.go(0)", tags$i("reset page content"), alt="Reset page"),
      tipify(downloadButton('downloadTable', 'Filtered Table'),"Download the filtered Statistical Table subset"),
      tipify(downloadButton('downloadIDList', 'Filtered ID-list'),"Download the filtered ENS geneID list for use in fpkm2heatmap or any gene list aware tool"),
      tipify(downloadButton('downloadMM', 'M&M'),"Download a text including the names and versions of all packages used in this webtool"),
      width=2
    ),

    # Show a filtering parameters and results
    mainPanel(
      textOutput('full.data.cnt'),
      textOutput("contrasts"),
      textOutput('min.lfc'),
      textOutput('max.pv'),
      textOutput('stringency'),
      br(),
      textOutput('filt.data.cnt'),
      DT::dataTableOutput("filt.data.table")
    )
  # end sidebarLayout
  )

  # end ui block
  )

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$downloadData <- downloadHandler(
    filename <- function() { paste("StatisticalResults", "xlsx", sep=".") },
    content <- function(file) { file.copy("Data/StatisticalResults.xlsx", file) }
  )

  load.data <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return("Waiting for data")

    # load data from excel file
    dat <- read.xlsx(inFile$datapath, sheet=1, rowNames=FALSE)

    # count contasts
    chr.col <- which(colnames(dat)==as.vector("Chromosome"))

    # store contrast names from LogFC column names
    contrast.names.all <- gsub(".:.logFC", "", colnames(dat)[seq(3, chr.col-1, 5)])

    # use the obtained values to populate the UI contrasts list
    output$choose_contrasts <- renderUI({
      # Create the checkboxes and select them all by default
      checkboxGroupInput("contrasts", "Check contrasts to be filtered",
                         choices  = contrast.names.all,
                         selected = contrast.names.all)
    })

    # reduce data size and keep only filterable columns
    data <- dat[,c(2, 1, sort(c(seq(3, chr.col-1, 5), seq(5, chr.col-1, 5))), chr.col:ncol(dat))]

    # return data as 'load.data()'
    data
  })

  output$full.data.cnt <- reactive({
    if (is.null(load.data())) return(NULL)
    paste("Rows in the Full data: ", nrow(load.data()))
  })

  output$contrasts <- renderText({
      if (is.null(load.data())) return(NULL)
      selected.contrasts <- paste(input$contrasts, collapse = ", ")
      paste("Selected contrasts (", length(input$contrasts), ") : ", selected.contrasts, sep="")
  })

  output$min.lfc <- renderText({
    paste("min log-FC (abs-val): ", input$min.lfc)
  })

  output$max.pv <- renderText({
    paste("max corrected-Pvalue: ", input$max.pv)
  })

  output$stringency <- renderText({
    paste("stringency: ", input$stringency)
  })

  output$keepunfiltered <- renderText({
    paste("Action for unfiltered columns is set to : ", input$keepunfiltered)
  })

  filtered.data <- eventReactive({input$goButton}, {
    # do nothing in absence of data
    if ( is.null(load.data()) | length(input$contrasts)==0 ) return(NULL)
    filtered <- load.data()

    # filter values in input$contrasts selected LR and pval columns
    #logfc.cols <- grepl( ".:.logFC" , names( filtered ) )
    #pval.cols <- grepl( ".:.FDR" , names( filtered ) )
    logfc.cols <- grepl( paste(paste(input$contrasts, ".:.logFC", sep=""), collapse="|") , names( filtered ) )
    pval.cols <- grepl( paste(paste(input$contrasts, ".:.FDR", sep=""), collapse="|") , names( filtered ) )

    # filter depending on input$stringency
    if (input$stringency=="one_or_more") {
      # handle case where only one contrast remains
      if (length(input$contrasts)>1) {
        filtered <- filtered[apply(abs(filtered[,logfc.cols]), 1, max) >= input$min.lfc,]
        filtered <- filtered[apply(filtered[,pval.cols], 1, min) <= input$max.pv,]
      } else {
        filtered <- filtered[abs(filtered[,logfc.cols]) >= input$min.lfc,]
        filtered <- filtered[filtered[,pval.cols] <= input$max.pv,]
      }
    } else {
      # input$stringency=="all_contrasts"
      # handle case where only one contrast remains
      if (length(input$contrasts)>1) {
        filtered <- filtered[apply(abs(filtered[,logfc.cols]), 1, min) >= input$min.lfc,]
        filtered <- filtered[apply(filtered[,pval.cols], 1, max) <= input$max.pv,]
      } else {
        filtered <- filtered[abs(filtered[,logfc.cols]) >= input$min.lfc,]
        filtered <- filtered[filtered[,pval.cols] <= input$max.pv,]
      }
    }

    # sort in increasing Pvalue order across all present contrasts
    # look at minimum pval in each row across all pval columns and sort rows based on these value
    # (Merci Rekin's)
    if (length(input$contrasts)>1) {
      filtered.sorted <- filtered[order(apply(filtered[,pval.cols], 1, min)),]
    } else {
      filtered.sorted <- filtered[order(filtered[,pval.cols]),]
    }

    # return final results
    filtered.sorted
    })

  output$filt.data.cnt <- reactive({
    if (is.null(filtered.data())) return(NULL)
    paste("gene rows in the filtered data: ", nrow(filtered.data()), sep="")
  })

  output$filt.data.table = DT::renderDataTable({
    if (is.null(filtered.data())) return(NULL)
    filtered.data()
  })

  output$downloadTable <- downloadHandler(
    filename = function() {
      paste(input$file1, "_filtered-data.xlsx", sep="")
    },
    content = function(file) {
      hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize=12,
                        fontName="Calibri", fgFill = "#4F80BD")
      write.xlsx(filtered.data(),
                 file,
                 sheetName = "Filtered-Data",
                 row.names = FALSE,
                 col.names = TRUE,
                 borders = "surrounding",
                 colWidths = "auto",
                 asTable = TRUE,
                 headerStyle = hs)
    }
  )

  output$downloadIDList <- downloadHandler(
    filename = function() {
      paste(input$file1, "_filteredID-list.txt", sep="")
    },
    content = function(file) {
      write(filtered.data()[,2], file = file, ncolumns = 1, sep = "\t")
    }
  )

  output$downloadMM <- downloadHandler(
    filename = function() {
      print(capture.output(paste(app.name, "_session_info.txt", sep="")))
    },
    content = function(file) {
      sink(file, append=TRUE)
      cat(paste("Thanks for using our tool", app.name, script.version, "\n", sep=" "))
      cat("\nYou can contact The Nucleomics Core at nucleomics@vib.be for any question\n")
      cat(paste("This data was generated on ", format(Sys.time(), "%a %b %d %H:%M:%S %Y"), "\n",sep=" "))
      cat("\nThe R packages used in the tools are listed next:\n")
      print(capture.output(sessionInfo()))
      sink()
    }
  )

  # end server block
}

# Run the application
shinyApp(ui = ui, server = server)

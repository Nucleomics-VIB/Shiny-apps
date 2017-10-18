# fpkm2heatmap.shinyapp
# A R/shiny tool to create a simple heatmap
# from a list of genes (signature)
# and a Nucleomics Core MS-Excel RNASeq count file

library("shiny")
library("shinyBS")
library("readr")
library("openxlsx")
library("DT")
library("pheatmap")
library("RColorBrewer")
library("grid")
library("ggplot2")
library("grDevices")

# you may uncomment the next line to allow large input files
# options(shiny.maxRequestSize=1000*1024^2)
# the following test checks if we are running on shinnyapps.io to limit file size dynamically
# ref: https://stackoverflow.com/questions/31423144/how-to-know-if-the-app-is-running-at-local-or-on-server-r-shiny/31425801#31425801
if ( Sys.getenv('SHINY_PORT') == "" ) { options(shiny.maxRequestSize=1000*1024^2) }

script.version="1.0"

# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application header
  headerPanel("Create a heamap plot for selected genes (RNASeq fpkm data)"),

  # Application title
  titlePanel(
    windowTitle = "RNASeq signature heatmap",
    tags$a(href="https://corefacilities.vib.be/nc", target="_blank",
           img(src='logo.png', align = "right", width="150", height="58.5", alt="VIB Nucleomics Core"))
  ),

  # Sidebar with input
  sidebarLayout(
    # show file import and molecule filters
    sidebarPanel(
      tags$h4(paste("code version: ", script.version, sep="")),
      tipify(downloadButton("downloadData", label = "Download test data"),
             "the Data is a MS-Excel file provided by the Nucleomics Core, with a fpkm worksheet reporting gene expression counts as second worksheet, you may produce a compatible file based on the test data provided here."),
      tipify(downloadButton("downloadSignature", label = "Download test signature"),
             "the signature is a two-lines text files with a first line starting by # followed by a space and a signature title (no spaces!), and a comma-separated list of Gene identifiers on line #2. These selected genes will be used to make a heatmap for all available samples."),
      tags$br(),
      tags$a(href="license.pdf", target="_blank", "usage licence"),
      tags$hr(),
      fileInput('file1', 'Choose RNASeq XLSX File', accept='.xlsx'),
      fileInput('file2', 'Choose text signature File', accept='.txt'),
      tags$h4("Edit settings & click ", tags$em("Plot")),
      textInput('outfile', "name for output File:", value="my_heatmap"),
      textInput('title', "Plot Title:", value="Custom HeatMap"),
      checkboxInput("log.trans", "Log2 transform (after adding 0.001)", value=TRUE),
      checkboxInput("show.gene.names", "Show Gene names:", value = TRUE),
      checkboxInput("show.sample.names", "Show Sample names:", value = TRUE),
      checkboxInput("show.legend", "Show legend:", value = TRUE),
      selectInput("drows", "Distance for genes:", c("NULL", "euclidean", "maximum", "manhattan"), selected="euclidean"),
      selectInput("dcols", "Distance for samples:", c("NULL", "euclidean", "maximum", "manhattan"), selected="euclidean"),
      selectInput("clustmet", "Clustering method:", c("average", "ward.D", "complete"), selected="average"),
      selectInput("color", "Color:", c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges",
                                       "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds",
                                       "YlGn", "YlGnBu", "YlOrBr", "YlOrRd")),
      selectInput("format", "Output format (png or pdf):", c("png", "pdf"), selected="png"),
      actionButton(inputId='goButton', "Plot", style='padding:4px; font-weight: bold; font-size:150%'),
      downloadButton('downloadTable', 'Download table'),
      downloadButton('downloadPlot', 'Download Plot')
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput('heatmap', width = "100%"),
      br(),
      textOutput('full.data.cnt'),
      textOutput('filt.data.cnt'),
      br(),
      DT::dataTableOutput("filt.data.table")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$downloadData <- downloadHandler(
    filename <- function() { paste("expXXXX-RNAseqCounts", "xlsx", sep=".") },
    content <- function(file) { file.copy("Data/expXXXX-RNAseqCounts.xlsx", file) }
  )

  output$downloadSignature <- downloadHandler(
    filename <- function() { paste("test.signature", "txt", sep=".") },
    content <- function(file) { file.copy("Data/test.signature.txt", file) }
  )

  fpkm.data <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)

    # load data from excel file
    dat <- read.xlsx(inFile$datapath, sheet=2)

    # keep only data columns (remove last columns including "Chromosome")
    chromosome.col <- which(colnames(dat)==as.vector("Chromosome"))
    fpkm.data <- dat[,c(2,1,3:(chromosome.col-1))]

    # remove some columns
    row.names(fpkm.data) <- paste(fpkm.data[,1], fpkm.data[,2], sep=":")
    fpkm.data <- fpkm.data[,-1]

    # kick useless part of names for samples
    colnames(fpkm.data) <- sub("@.*", "", colnames(fpkm.data))

    # return data as 'fpkm.data()'
    fpkm.data
  })

  output$full.data.cnt <- reactive({
    if (is.null(fpkm.data())) return(NULL)
    paste("Rows in the Full data: ", nrow(fpkm.data()))
  })

  sig.name <- reactive({
    inFile <- input$file2
    if (is.null(inFile)) return(NULL)

    # read signature name
    sig.name <- readLines(inFile$datapath, n=1)

    # check format or die
    if( startsWith(sig.name, "# ") ) {
      sig.name <- gsub("# ", "", sig.name)
    } else {
      stop("Signature first row should be '# signature_name' (only space between # and name)")
    }
    sig.name
  })

  sig.vect <- reactive({
    inFile <- input$file2
    if (is.null(inFile)) return(NULL)

    # read signature ID list
    sig.vect <- read.table(inFile$datapath, skip=1, sep=",", header=FALSE)
    sig.vect <- as.vector(t(sig.vect))

    if(length(sig.vect)==0) {
      stop("Signature second row should be a comma-separated list of ENSEMBL-IDs")
    }

    sig.vect
    })

  filtered.data <- eventReactive(input$goButton, {
    # do nothing in absence of data
    if (is.null(fpkm.data())) return(NULL)
    if (is.null(sig.vect())) return(NULL)

    # select only signature rows and discard Gene.ID column to keep only FPKM in data.frame
    fpkm.data <- as.data.frame(fpkm.data())
    hm.data <- fpkm.data[fpkm.data$Gene.ID %in% sig.vect(), 2:length(fpkm.data)]

    # return data
    hm.data
    })

  output$filt.data.cnt <- reactive({
    paste("Rows in the Signature data: ", nrow(filtered.data()))
  })

  output$filt.data.table = DT::renderDataTable({
    if (is.null(filtered.data())) return(NULL)

    hm.data <- filtered.data()

    if (input$log.trans==TRUE) {
      hm.data <- round(log(hm.data+0.001, 2),3)
    }

    hm.data

  })

  plotInput <- function(){
    if (is.null(filtered.data())) return(NULL)

    # define metrics for clustering
    if (input$drows=="NULL") {
      cluster.rows=FALSE
      drows=NULL
    } else {
      cluster.rows=TRUE
      drows=input$drows
    }

    if (input$dcols=="NULL") {
      cluster.cols=FALSE
      dcols=NULL
    } else {
      cluster.cols=TRUE
      dcols=input$dcols
    }

    if (input$clustmet=="NULL") {
      clustmet=NULL
    } else {
      clustmet=input$clustmet
    }

    hm.data <- filtered.data()

    if (input$log.trans==TRUE) {
      hm.data <- round(log(hm.data+0.001, 2),3)
    }

    heatmap <- pheatmap(hm.data,
                        color = (brewer.pal(9, input$color)),
                        fontsize = 8,
                        cellwidth = 12, cellheight = 12, scale = "none",
                        treeheight_row = 100,
                        kmeans_k = NA,
                        show_rownames = input$show.gene.names,
                        show_colnames = input$show.sample.names,
                        main = input$title,
                        clustering_method = clustmet,
                        cluster_rows = cluster.rows,
                        cluster_cols = cluster.cols,
                        clustering_distance_rows = drows,
                        clustering_distance_cols = dcols,
                        legend = input$show.legend)
    heatmap
  }

  output$heatmap <- renderPlot({
    print(plotInput())
  })

  output$downloadPlot <- downloadHandler(
    filename =  function() {
      paste(input$outfile, input$format, sep=".")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      if(input$format == "png")
        png(file) # open the png device
      else
        pdf(file, onefile=FALSE) # open the pdf device
      # plot
      plotInput()
      dev.off()  # turn the device off
    })

  output$downloadTable <- downloadHandler(
    filename = function() {
      paste(input$outfile, ".xlsx", sep="")
    },
    content = function(file) {
      write.xlsx(filtered.data(), file, row.names=TRUE, col.names=TRUE)
    }
  )

}

# Run the application
shinyApp(ui = ui, server = server)

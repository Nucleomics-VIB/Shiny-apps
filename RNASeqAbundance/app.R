# rnaseqabundance.shinyapp
# A R/shiny tool to plot RNASeq abundance plots
# from a Nucleomics Core MS-Excel RNASeq count file

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

app.name <- "RNASeqAbundance"
script.version <- "1.0.0"

# maximum signature length
maxlen <- 200

# Define UI for application that draws a histogram
ui <- fluidPage(
  HTML('<style type="text/css">
    .row-fluid { width: 20%; }
       .well { background-color: #99CCFF; }
       .shiny-html-output { font-size: 14px; line-height: 15px; }
       </style>'),
  # Application header
  headerPanel("Create a plot of gene abundance from RNASeq Raw Counts data"),

  # Application title
  titlePanel(
    windowTitle = "RNASeq Abundance plot",
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
      tipify(fileInput('file1', 'Choose RNASeq XLSX File', accept='.xlsx'),
             "the Data is a MS-Excel file provided by the Nucleomics Core, with worksheet#1 reporting Raw Counts, you may produce a compatible file based on the test data provided here."),
      tags$h4("Edit settings & click ", tags$em("Plot")),
      actionButton(inputId='goButton', "Plot", style='padding:4px; font-weight: bold; font-size:150%'),
      tipify(selectInput("genename", "Gene name:", c("Gene_symbol", "ENSembl_GID", "both"), selected="Gene_symbol"),"the gene names to show at the right of the rows."),
      textInput('title', "Plot Title:", value="Gene Abundance Plot"),
      selectInput("format", "Output format (png or pdf):", c("png", "pdf"), selected="png"),
      textInput('outfile', "name for output File:", value="my_plot"),
      downloadButton('downloadTable', 'Download table'),
      downloadButton('downloadPlot', 'Download Plot'),
      tipify(downloadButton('downloadMM', 'Download M&M'),"Download a text including the names and versions of all packages used in this webtool")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput('heatmap', width = "100%"),
      br(),
      textOutput('full.data.cnt'),
      br(),
      div(DT::dataTableOutput("full.data.table"), style = "font-size: 75%; width: 75%")
    )
  )

  # end UI block
  )

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$downloadData <- downloadHandler(
    filename <- function() { paste("expXXXX-RNAseqCounts", "xlsx", sep=".") },
    content <- function(file) { file.copy("Data/expXXXX-RNAseqCounts.xlsx", file) }
  )

  count.data <- reactive({
    inFile <- input$file1

    if (is.null(inFile)) return(NULL)

    # load data from excel file
    dat <- read.xlsx(inFile$datapath, sheet=1)

    # keep only data columns (remove last columns including "Chromosome")
    chromosome.col <- which(colnames(dat)==as.vector("Chromosome"))
    count.data <- dat[,c(2,1,3:(chromosome.col-1))]

	# define row names
    # row.names(fpkm.data) <- paste(fpkm.data[,1], fpkm.data[,2], sep=":")
	# define gene names based on input$genename ("Gene_symbol", "ENSembl_GID", "both")
	if ( input$genename == "Gene_symbol") {
		row.names(count.data) <- count.data[,1]
		} else if ( input$genename == "ENSembl_GID") {
			row.names(count.data) <- count.data[,2]
			} else { 
				row.names(count.data) <- paste(count.data[,1], count.data[,2], sep=":")
				}
    # remove some columns
    count.data <- count.data[,-1]

    # kick useless part of names for samples
    colnames(count.data) <- sub("@.*", "", colnames(count.data))

    # return data as 'fpkm.data()'
    count.data
  })

  output$full.data.cnt <- reactive({
    if (is.null(count.data())) return("Waiting for data!")

    paste("Rows in the Full data: ", nrow(count.data()))
  })

  output$full.data.table = DT::renderDataTable({
    if (is.null(count.data())) return(NULL)

    count.data()
  })

  
  output$downloadPlot <- downloadHandler(
    filename =  function() {
      paste(input$outfile, input$format, sep=".")
    },
    # content is a function with argument file. content writes the plot to file
    content = function(file) {
      do.call("pheatmap", c(hm.parameters(), filename=file))
    })

  output$downloadTable <- downloadHandler(
    filename = function() { paste(input$outfile, "_data.xlsx", sep="") },
    content = function(file) {
      hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize=12,
                        fontName="Calibri", fgFill = "#4F80BD")
      write.xlsx(filtered.data(),
                 file,
                 sheetName = "Heatmap_Data",
                 row.names = TRUE,
                 col.names = TRUE,
                 borders = "surrounding",
                 colWidths = "auto",
                 asTable = TRUE,
                 headerStyle = hs)
    }
  )

  output$downloadMM <- downloadHandler(
    filename = function() {
      paste(input$outfile, "_session_info.txt", sep="")
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

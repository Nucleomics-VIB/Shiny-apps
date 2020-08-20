# fpkm2heatmap.shinyapp
# A R/shiny tool to create a simple heatmap
# from a list of genes (signature)
# and a Nucleomics Core MS-Excel RNASeq count file

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

app.name <- "fpkm2heatmap"
script.version <- "1.4.0"

# maximum signature length
maxlen <- 500
deflen <- 50

# Define UI for application that draws a histogram
ui <- fluidPage(
  HTML('<style type="text/css">
    .row-fluid { width: 20%; }
       .well { background-color: #99CCFF; }
       .shiny-html-output { font-size: 14px; line-height: 15px; }
       </style>'),
  # Application header
  headerPanel("Create a heatmap plot for selected genes (RNASeq fpkm data)"),

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
      tags$h5(paste(app.name, " version: ", script.version, sep="")),
      downloadButton("downloadData", label = "Download test data"),
      downloadButton("downloadSignature", label = "Download test signature"),
      tags$br(),
      tags$a(href="license.pdf", target="_blank", "usage licence"),
      tags$br(),
      tags$a(href="javascript:history.go(0)", tags$i("reset page content"), alt="Reset page"),
      tags$hr(),
      tipify(fileInput('file1', 'Choose RNASeq XLSX File', accept='.xlsx'),
             "the Data is a MS-Excel file provided by the Nucleomics Core, with worksheet#2 reporting gene expression (FPKM), you may produce a compatible file based on the test data provided here."),
      tipify(fileInput('file2', 'Choose text signature File', accept='.txt'),
             "the signature is a one-column text files containing one EnsEMBL ID per line. It can for instance be made from the top-N DE genes in your data or from a list of genes members of a pathway or biological function. We limit here the length of a signature to 200 to prevent generating plots that cannot be printed on one page (the first 200 IDs are used if the list is larger)"),
      tags$h4("Edit settings & click ", tags$em("Plot")),
      actionButton(inputId='goButton', "Plot", style='padding:4px; font-weight: bold; font-size:150%'),
      tipify(selectInput("genename", "Gene name:", c("Gene_symbol", "ENSembl_GID", "both"), selected="Gene_symbol"),"the gene names to show at the right of the rows."),
      textInput('title', "Plot Title:", value="Custom HeatMap"),
      sliderInput("obs", "Number of genes to plot: ", min=1, max=maxlen, value=deflen),
      checkboxInput("show.gene.names", "Show Gene names:", value = TRUE),
      checkboxInput("show.sample.names", "Show Sample names:", value = TRUE),
      checkboxInput("show.legend", "Show legend:", value = TRUE),
      tipify(selectInput("scale.data", "Scale data:", c("none", "row", "column"), selected="row"),"the values should be centered and scaled in either the row direction or the column direction, or none."),
      tipify(selectInput("drows", "Distance for genes:", c("none", "euclidean", "maximum", "manhattan", "camberra", "binary", "minkowski"), selected="euclidean"),"distance measure to compute the distances between the genes of the fpkm matrix"),
      tipify(selectInput("dcols", "Distance for samples:", c("none", "euclidean", "maximum", "manhattan", "camberra", "binary", "minkowski"), selected="euclidean"),"distance measure to compute the distances between the samples of the fpkm matrix"),
      tipify(selectInput("clustmet", "Clustering method:", c("average", "ward.D", "complete","ward.D2","single"), selected="average"),"the agglomeration method to be used"),
      selectInput("color", "Color:", c("BlueWhiteRed","Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges",
                                       "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds",
                                       "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"), selected="RedWhiteBlue"),
      selectInput("format", "Output format (png or pdf):", c("png", "pdf"), selected="png"),
      textInput('outfile', "name for output File:", value="my_heatmap"),
      downloadButton('downloadTable', 'Download table'),
      downloadButton('downloadPlot', 'Download Plot'),
      tipify(downloadButton('downloadMM', 'Download M&M'),"Download a text including the names and versions of all packages used in this webtool")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput('heatmap', width = "100%"),
      br(),
      textOutput('full.data.cnt'),
      textOutput('filt.data.cnt'),
      br(),
      div(DT::dataTableOutput("filt.data.table"), style = "font-size: 75%; width: 75%")
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

	# define row names
    # row.names(fpkm.data) <- paste(fpkm.data[,1], fpkm.data[,2], sep=":")
	# define gene names based on input$genename ("Gene_symbol", "ENSembl_GID", "both")
	if ( input$genename == "Gene_symbol") {
		row.names(fpkm.data) <- make.unique(fpkm.data[,1],sep = ".")
		} else if ( input$genename == "ENSembl_GID") {
			row.names(fpkm.data) <- fpkm.data[,2]
			} else { 
				row.names(fpkm.data) <- paste(fpkm.data[,1], fpkm.data[,2], sep=":")
				}
    # remove some columns
    fpkm.data <- fpkm.data[,-1]

    # kick useless part of names for samples
    colnames(fpkm.data) <- sub("@.*", "", colnames(fpkm.data))

    # return data as 'fpkm.data()'
    fpkm.data
  })

  output$full.data.cnt <- reactive({
    if (is.null(fpkm.data())) return("Waiting for data!")

    paste("Rows in the Full data: ", nrow(fpkm.data()))
  })

  sig.vect <- reactive({
    inFile <- input$file2

    if (is.null(inFile)) return(NULL)

        # read signature ID list
    sig.vect <- read.table(inFile$datapath, sep=",", header=FALSE)
    sig.vect <- as.vector(t(sig.vect))

    if(length(sig.vect)==0) {
      stop("Signature should be a 1-column text file with EnsEMBL IDs")
    }

    sig.vect
    })

  filtered.data <- eventReactive({input$goButton | input$obs}, {
    # do nothing in absence of data
    if (is.null(fpkm.data())) return(NULL)
    if (is.null(sig.vect())) return(NULL)

    # select only signature rows and discard Gene.ID column to keep only FPKM in data.frame

    # limit to maxlen
    sig.vect <- sig.vect()[1:input$obs]

    fpkm.data <- as.data.frame(fpkm.data())
    hm.data <- fpkm.data[fpkm.data$Gene.ID %in% sig.vect, 2:length(fpkm.data)]

    # return data
    hm.data
    })

  output$filt.data.cnt <- reactive({
    if (is.null(filtered.data())) return("Waiting for data!")
    paste("Rows in the Signature data: ", nrow(filtered.data()))
  })

  output$filt.data.table = DT::renderDataTable({
    if (is.null(filtered.data())) return(NULL)

    filtered.data()
  })

  hm.parameters <- function(){
    if (is.null(filtered.data())) return(NULL)

    # define metrics for clustering
    if (input$drows=="none") {
      cluster.rows=FALSE
      drows=NULL
    } else {
      cluster.rows=TRUE
      drows=input$drows
    }

    if (input$dcols=="none") {
      cluster.cols=FALSE
      dcols=NULL
    } else {
      cluster.cols=TRUE
      dcols=input$dcols
    }

    # only active when at least one above is set
    clustmet=input$clustmet

    hm.data <- filtered.data()

    # about colors
    if (input$color=="BlueWhiteRed") {
      sel.pal <- colorRampPalette(c("blue","white","red"))(256)
    } else {
      sel.pal <- brewer.pal(9, input$color)
    }

    hm.parameters <- list(hm.data,
                        color = sel.pal,
                        fontsize = 8,
                        cellwidth = 12, cellheight = 12,
                        treeheight_row = 100,
                        kmeans_k = NA,
                        show_rownames = input$show.gene.names,
                        show_colnames = input$show.sample.names,
                        main = input$title,
                        scale = input$scale.data,
                        clustering_method = clustmet,
                        cluster_rows = cluster.rows,
                        cluster_cols = cluster.cols,
                        clustering_distance_rows = drows,
                        clustering_distance_cols = dcols,
                        legend = input$show.legend)
    hm.parameters
  }

  output$heatmap <- renderPlot({
    if (is.null(hm.parameters())) return(NULL)
    do.call("pheatmap", hm.parameters())
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
      cat ("This tool generates a Heatmap after clustering the gene expression values (FPKM), or both FPKM values and samples.")
      cat (" Different distance metrics and options can be set by the user although defaults should do for most cases.")
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

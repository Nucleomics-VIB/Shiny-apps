# assemblyNplot.shinyapp
# a Shiny web application that draws N graphs from a zip of denovo assemblies
library("shiny")
library("shinyBS")
library("seqinr")
library("ggplot2")
library("scales")
library("DT")

# you may uncomment the next line to allow large input files (4GB)
options(shiny.maxRequestSize=4*1000*1024^2)
# the following test checks if we are running on shinnyapps.io to limit file size dynamically
# ref: https://stackoverflow.com/questions/31423144/how-to-know-if-the-app-is-running-at-local-or-on-server-r-shiny/31425801#31425801
#if ( Sys.getenv('SHINY_PORT') == "" ) { options(shiny.maxRequestSize=1000*1024^2) }

app.name <- "assemblyNplot"
script.version <- "1.1"

# cleanup from previous uploads
cleanup <- function () {
  folders <- list.dirs('.', recursive=FALSE)
  # keep only the following folders
  keep <- c("www", "Data")
  remove <- subset(folders, !grepl(paste0(keep, collapse="|"), folders))
  unlink(remove, recursive=TRUE)
}

# measure sequence lengths from fasta
Fasta2length <-function(fastaFile) {
  #fa <- read.fasta(file = fastaFile, as.string = TRUE, seqonly = TRUE)
  fl <- getLength(read.fasta(file = fastaFile, as.string = TRUE, seqonly = TRUE))
  # return a vector of lengths
  as.vector(fl)
}

# compute XNXX function
Nvalue <- function(lim, x, na.rm = TRUE) {
  # handle NA values
  if(isTRUE(na.rm)){
    x <- x[!is.na(x)]
    }
  cutval <- 100/lim
  # compute LXX and NXX
  sorted <- sort(x, decreasing = TRUE)
  SXX <- sum(x)/cutval
  csum <- cumsum(sorted)
  GTLXX <- as.vector(csum >= SXX)
  LXX=min(which(GTLXX == TRUE))
  NXX <- round(sorted[LXX], 1)
  # eg: get NXX with lst['NXX']
  NXX
  }

# format with thousand separator
fnum <- function(x) {
  return(format(as.numeric(x), nsmall=0, big.mark="'"))
}

# CLEANUP OLD DATA
cleanup()
    
# Define UI for application that draws a histogram
ui <- fluidPage(
  HTML('<style type="text/css">
       .row-fluid { width: 25%; }
       .well { background-color: #99CCFF; }
       .shiny-html-output { font-size: 14px; line-height: 15px; }
       </style>'),
  
  # Application header
  headerPanel("Plot N-graphs from a Zip of de-novo assembly fasta files"),
  
  # Application title
  titlePanel(
    windowTitle = "AssemblyNplot",
    tags$a(href="https://corefacilities.vib.be/nc", target="_blank",
           img(src='logo.png', align = "right", width="150", height="58.5", alt="VIB Nucleomics Core"))
    ),
  
  sidebarLayout(
    # show file import and molecule filters
    sidebarPanel(
      tags$h5(paste(app.name, " version: ", script.version, sep="")),
      downloadButton("downloadData", label = "Download test data"),
      tags$br(),
       tipify(fileInput("upload", "Upload", accept = ".zip"), 
              "A zip files containing all fasta assemblies to plot"),
       br(),
       actionButton("process", "Process uploaded data"),
       hr(),
       textInput('outfile', "name for output File:", value="assemblyNplot"),
       selectInput("format", "Output format (png or pdf):", c("png", "pdf"), selected="pdf"),
       downloadButton('downloadPlot', 'Download Plot')
    ),
  
  mainPanel(
    plotOutput('plot', width = "100%"),
    div(DT::dataTableOutput('ntable'), style = "font-size: 75%; width: 75%")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$downloadData <- downloadHandler(
    filename <- function() { "testData.zip" },
    content <- function(file) { file.copy("www/testData.zip", file) }
  )
  
  fasta.files <- eventReactive({input$process}, {
    # remove previous uploads
    unlink("Data", recursive=TRUE)
    # unzip user data
    unzip.files <- unzip(input$upload$datapath, list = FALSE, exdir = "Data")
    # get rid of OSX hidden and empty stuff
    fasta.files <- subset(unzip.files, !grepl("__MACOSX|.DS_Store|/$", unzip.files))
    unlink("Data/__MACOSX", recursive=TRUE)
    # return fasta fiel list
    fasta.files
    })
  
  parse.data <- reactive({
    if (is.null(fasta.files())) return(NULL)
    
    # initialize
    n.table <- data.frame()
    n <- length(fasta.files())
    
    withProgress(message = 'Analyzing ', value = 0, {
      for (assembly in fasta.files()){
        title <- basename(assembly)
        incProgress(1/n, detail = title)
        lengths <- Fasta2length(assembly)
        width <- sum(lengths)
        x <- seq(1, 100, by=1)
        y <- sapply(x, function(x) Nvalue(x, lengths))
        name <- rep(title, length(x))
        width <- rep(fnum(width), length(x))
        dat <- data.frame(assembly=paste0(name," (", width,")"), x=x, y=y)
        n.table <- rbind(n.table, dat)
        }
      })
    as.data.frame(n.table)
  })

  output$ntable = DT::renderDataTable({
    if (is.null(parse.data())) return(NULL)
    parse.data()
    })
  
  plotInput <- reactive({
    df <- parse.data()
    # convert to kilobases in plot
    p <- ggplot(data=df, aes(x=x, y=y/1000, group=assembly, colour=assembly)) + 
      scale_y_continuous(trans="log10",labels = waiver()) +
      annotation_logticks(sides = "l") +
      geom_line(size = 0.75, linetype="dotted") + 
      geom_point(aes(shape=assembly), size = 2) +
      geom_vline(xintercept = 50, linetype="dotted", 
                 color = "red", size=0.5) +
      ggtitle("NG graphs of the assemblies in scaffold length") + 
      labs(x = "NG%", y = "Contig/Scaffold length (kb)") +
      labs(x = "NG", y = "Scaffold NG length (kb)") +
      theme(axis.text.x = element_text(colour="grey20",size=8,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=8,angle=0,hjust=1,vjust=0,face="plain"),
            axis.title.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=0,face="plain"),
            axis.title.y = element_text(colour="grey20",size=10,angle=90,hjust=.5,vjust=.5,face="plain"),
            legend.justification = c(0,1),
            legend.position = c(0.1,0.5),
            legend.title = element_blank(),
            legend.text = element_text(size=10),
            legend.key = element_rect(colour = NA, fill = NA),
            legend.key.size = unit(0.8, 'lines'),
            legend.background = element_rect(fill="transparent"),
            plot.title = element_text(margin=margin(b=0), size = 16))
  })
  
  output$plot <- renderPlot({
    plot(plotInput(), width="640px", height="480px")
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste(input$outfile, input$format, sep=".") },
    content = function(file) {
      if(input$format == "png")
        png(file, width = 640, height = 480, units = "px") # open the png device
      else
        pdf(file, width = 8, height = 6) # open the pdf device
      plot(plotInput())
      dev.off()  # turn the device off
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)

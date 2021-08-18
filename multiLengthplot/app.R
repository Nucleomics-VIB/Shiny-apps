# multiLengthPlot.shinyapp
# a Shiny web application that draws density plot from a zZip archive of length distributions
library("shiny")
library("shinyBS")
library("ggplot2")
library("scales")
library("DT")

# you may uncomment the next line to allow large input files (4GB)
options(shiny.maxRequestSize=1000*1024^2)
# the following test checks if we are running on shinnyapps.io to limit file size dynamically
# ref: https://stackoverflow.com/questions/31423144/how-to-know-if-the-app-is-running-at-local-or-on-server-r-shiny/31425801#31425801
#if ( Sys.getenv('SHINY_PORT') == "" ) { options(shiny.maxRequestSize=1000*1024^2) }

app.name <- "multiLengthPlot"
script.version <- "1.4"

# cleanup from previous uploads
cleanup <- function () {
  folders <- list.dirs('.', recursive=FALSE)
  # keep only the following folders
  keep <- c("www", "Data")
  remove <- subset(folders, !grepl(paste0(keep, collapse="|"), folders))
  unlink(remove, recursive=TRUE)
}

# compute XNXX function
Nvalue <- function(lim, x, na.rm = TRUE){
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
fnum <- function(x){
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
  headerPanel(paste0(
    "Plot from a Zip of length distributions (v",
    script.version,")")),
  
  # Application title
  titlePanel(
    windowTitle = "multiLengthPlot",
    tags$a(href="https://corefacilities.vib.be/nc", target="_blank",
           img(src='logo.png', align = "right", width="150", height="58.5", alt="VIB Nucleomics Core"))
  ),
  
  sidebarLayout(
    sidebarPanel(
      tags$h5(paste(app.name, " version: ", script.version, sep="")),
      # show file import and length filters
      downloadButton("downloadData", label = "Download test data"),
      tags$br(),
      tipify(fileInput('upload', 'Upload', accept = c('.zip')), 
             "A zip files containing 1 or more length distributions to plot"),
      tipify(textInput('maxrec', 'max records per file', value=100000),
      "Do not exceed 5M to avoid reaching the RAM limit"),
      radioButtons('xscale', 'X-Scale',
                   choices = c(Linear = "lin",
                               Log = "log"),
                   selected = "lin"),
      radioButtons('stat', 'Stat',
                   choices = c(density = "density",
                               scaled = "scaled"),
                   selected = "density"),
      textInput('minl', 'Min length', value=0),
      textInput('maxl', 'Max length', value=1E+5),
      actionButton('process', 'Plot filtered data'),
      hr(),
      textInput('outfile', 'name for output File:', value="densityPlot"),
      selectInput('format', 'Output format (png or pdf):', c("png", "pdf"), selected="pdf"),
      downloadButton('downloadPlot', 'Download Plot')
    ),
    
    mainPanel(
      plotOutput('plot', width = "100%"),
      div(DT::dataTableOutput('infotable'), style = "font-size: 75%; width: 75%; align: right")
    )
  )
  )

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$downloadData <- downloadHandler(
    filename <- function() { "testData.zip" },
    content <- function(file) { file.copy("www/testData.zip", file) }
  )
  
  density.files <- eventReactive({input$process}, {
    # remove previous uploads
    unlink("Data", recursive=TRUE)
    # unzip user data
    unzip.files <- unzip(input$upload$datapath, list = FALSE, exdir = "Data")
    # get rid of OSX hidden and empty stuff
    density.files <- subset(unzip.files, !grepl("__MACOSX|.DS_Store|/$", unzip.files))
    unlink("Data/__MACOSX", recursive=TRUE)
    # return file list
    density.files
  })
  
  import.data <- reactive({
    if (is.null(density.files())) return(NULL)
    # initialize
    n <- length(density.files())
    data <- data.frame()
    withProgress(message = 'Importing ', value = 0, {
      for (dfile in density.files()){
        title <- basename(dfile)
        title <- gsub('.txt$','',title)
        incProgress(1/n, detail = title)
        # collect lengths from a single file with scan
        all.lengths <- scan(dfile, 
                        numeric(), 
                        quote = "", 
                        blank.lines.skip = TRUE)
        # data size
        recnum <- length(all.lengths)

        # take a random sample of smallest(maxrec, recnum)
        sample.size <- min(as.numeric(c(recnum, input$maxrec)))
        lengths <- sample(all.lengths, sample.size)
        
        dat <- data.frame(name=title, len=as.vector(lengths))
        data <- rbind(data, dat)
      }
    })
    return(data)
  })
  
  filter.data <- reactive({
    if (is.null(import.data())) return(NULL)
    
    # filter data based on limits
    data <- import.data()
    minl <- as.numeric(input$minl)
    maxl <- as.numeric(input$maxl)    
    data <- subset(data, len>minl & len<maxl)
    
    # get info for each dataset
    datasets <- as.vector(unique(data$name))
    n <- length(datasets)
    info <- data.frame()
    withProgress(message = 'Analyzing ', value = 0, {
      for (ds in datasets){
        incProgress(1/n, detail = ds)
        # collect lengths from a single file
        ds.data <- subset(data, data$name==ds)
        # collect metrics
        nrec <- nrow(ds.data)
        width <- sum(ds.data$len)
        n50 <- Nvalue(50, ds.data$len)
        dat <- data.frame(name=ds, nrec=nrec, width=width, n50=n50)
        info <- rbind(info, dat)
      }
    })
    return(list(data=data, info=info))
  })
  
  output$infotable = DT::renderDataTable({
    if (is.null(filter.data())) return(NULL)
    info <- filter.data()$info
    info[,2:4] <- format(info[,2:4],nsmall = 0, big.mark = "'")
    info
  }, option=list(dom = 't', 
                 columnDefs=list(list(targets=2:4, class="dt-right"))
                 )
  )
  
  plotInput <- reactive({
    if (is.null(filter.data())) return(NULL)
    info <- filter.data()$info
    df <- filter.data()$data
    
    if(input$stat == "scaled") {
      p <- ggplot(data=df, aes(x=len, y=..scaled.., group=name, colour=name)) +
        ggtitle("Frequency distributions") + 
        labs(x = "length", y = "frequency")
    } else {
      p <- ggplot(data=df, aes(x=len, group=name, colour=name)) +
        ggtitle("Density distributions") + 
        labs(x = "length", y = "density")
    }
    
    p <- p + geom_density(size=1) + 
      theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=10,angle=0,hjust=1,vjust=0,face="plain"),
            axis.title.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=0,face="plain"),
            axis.title.y = element_text(colour="grey20",size=10,angle=90,hjust=.5,vjust=.5,face="plain"),
            legend.justification = c(0,1),
            legend.title = element_blank(),
            legend.text = element_text(size=12),
            legend.key = element_rect(colour = NA, fill = NA),
            legend.key.size = unit(0.8, 'lines'),
            legend.background = element_rect(fill="transparent"),
            plot.title = element_text(margin=margin(b=0), size = 14))
    
    if(input$xscale == "log") {
      p <- p + scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      ) + annotation_logticks(sides="b")
    }
    
    # add N50 lines and plot
    p + geom_vline(data=info, aes(xintercept=n50, group=name, colour=name), size=0.75)
  })
  
  output$plot <- renderPlot({
    plot(plotInput(), width="640px", height="480px")
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste(input$outfile, input$format, sep=".") },
    content = function(file) {
      if(input$format == "png")
        png(file, width = 800, height = 480, units = "px") # open the png device
      else
        pdf(file, width = 10, height = 6) # open the pdf device
      plot(plotInput())
      dev.off()  # turn the device off
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
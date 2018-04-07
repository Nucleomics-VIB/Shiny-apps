# LengthPlot.shinyapp
# a Shiny web application that draws a length distribution plot from a text file of lengths

library("shiny")
library("shinyBS")
library("ggplot2")
library("scales")

# you may uncomment the next line to allow large input files
options(shiny.maxRequestSize=1000*1024^2)
# the following test checks if we are running on shinnyapps.io to limit file size dynamically
# ref: https://stackoverflow.com/questions/31423144/how-to-know-if-the-app-is-running-at-local-or-on-server-r-shiny/31425801#31425801
#if ( Sys.getenv('SHINY_PORT') == "" ) { options(shiny.maxRequestSize=1000*1024^2) }

app.name <- "LengthPlot"
script.version <- "1.0"

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

# Define UI for application that draws a histogram
ui <- fluidPage(
  HTML('<style type="text/css">
       .row-fluid { width: 25%; }
       .well { background-color: #99CCFF; }
       .shiny-html-output { font-size: 14px; line-height: 15px; }
       </style>'),
  
  # Application header
  headerPanel("Plot the length distribution from a text file"),
  
  # Application title
  titlePanel(
    windowTitle = "LengthPlot",
    tags$a(href="https://corefacilities.vib.be/nc", target="_blank",
           img(src='logo.png', align = "right", width="150", height="58.5", alt="VIB Nucleomics Core"))
  ),

  sidebarLayout(
    sidebarPanel(
      tipify(fileInput("upload", "Upload", accept = c('.txt', '.text')), 
             "A text file reporting lengths to plot"),
      radioButtons("xscale", "X-Scale",
                   choices = c(Linear = "lin",
                               Log = "log"),
                   selected = "lin"),
      textInput('minl', 'Min length', value=0),
      textInput('maxl', 'Max length', value=1E+5),
      actionButton("process", "Plot filtered data"),
      hr(),
      textInput('outfile', "name for output File:", value="LengthPlot"),
      selectInput("format", "Output format (png or pdf):", c("png", "pdf"), selected="pdf"),
      downloadButton('downloadPlot', 'Download Plot')
    ),
    
    mainPanel(
      textOutput("selected_var"),
      plotOutput('plot', width = "100%")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$selected_var <- renderText({ 
    paste("You have selected as limits: ", input$minl, " - ", input$maxl)
  })

  filter.data <- eventReactive({input$process}, {
    req(input$upload)
    # "filtered_albacore_reads_lengths.txt"
    # input$upload$datapath
    df <- read.table(input$upload$datapath, 
                   stringsAsFactors = FALSE,
                   header = FALSE,
                   sep = "",
                   col.names=c("len"),
                   colClasses=c("numeric"),
                   blank.lines.skip=TRUE
                   )
    df <- as.data.frame(df)
    minl <- as.numeric(input$minl)
    maxl <- as.numeric(input$maxl)    
    df <- subset(df, len>minl & len<maxl)
    df
  })
  
  plotInput <- reactive({
    if (is.null(filter.data())) return(NULL)
    df <- filter.data()
    width <- sum(df$len)
    N50 <- Nvalue(50, df$len)
    nseq <- nrow(df)
    # plot
    p <- ggplot(data=df, aes(x=len)) + 
      geom_density(col="blue", size=1) + 
      geom_vline(xintercept = N50, linetype="dotted", 
                 color = "red", size=1) +
      ggtitle(paste0("Sequence length distribution for ", fnum(nseq), " sequences (", fnum(width), " bases; N50=", fnum(N50), ")")) + 
      labs(x = "sequence length (bases)", y= "density") +
      theme(axis.text.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=10,angle=0,hjust=1,vjust=0,face="plain"),
            axis.title.x = element_text(colour="grey20",size=10,angle=0,hjust=.5,vjust=0,face="plain"),
            axis.title.y = element_text(colour="grey20",size=10,angle=90,hjust=.5,vjust=.5,face="plain"),
            legend.justification = c(0,1),
            legend.position = c(0.1,0.5),
            legend.title = element_blank(),
            legend.text = element_text(size=10),
            legend.key = element_rect(colour = NA, fill = NA),
            legend.key.size = unit(0.8, 'lines'),
            legend.background = element_rect(fill="transparent"),
            plot.title = element_text(margin=margin(b=0), size = 14))
    
    if(input$xscale == "log") {
      p <- p + scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        annotation_logticks(sides="b")
    }
    # return plot object
    p
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

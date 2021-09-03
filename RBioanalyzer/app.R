# RBioanalyzer.shinyapp
# A R/shiny tool to plot from Bioanalyzer exported text files
#
# Stephane Plaisance, VIB Nucleomics Core
# visit our Git: https://github.com/Nucleomics-VIB 
# version: 2017-03-17_v1.0
# © by using this tool, you accept the licence saved under ./www/licence.pdf

library("shiny")
library("readr")
library("stringr")
library("data.table")
library("lattice")
library("latticeExtra")

# you may uncomment the next line to allow large input files
options(shiny.maxRequestSize=1000*1024^2)
# the following test checks if we are running on shinnyapps.io to limit file size dynamically
# ref: https://stackoverflow.com/questions/31423144/how-to-know-if-the-app-is-running-at-local-or-on-server-r-shiny/31425801#31425801
#if ( Sys.getenv('SHINY_PORT') == "" ) { options(shiny.maxRequestSize=1000*1024^2) }

script.version="1.1 (2021-09-03)"

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application header 
  headerPanel(
    "R-Bioanalyzer overlay plot"
  ),
  
  # Application title
  titlePanel(
    windowTitle = "R-Bioanalyzer overlay plot",
    tags$a(href="https://corefacilities.vib.be/nc", target="_blank", 
           img(src='logo.png', align = "right", width="150", height="58.5", alt="VIB Nucleomics Core"))
  ),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      tags$p("Load two or three Bioanalyzer data files (csv) and create an overlay plot."),
      tags$p(tags$em("Tip: You may switch the ladder with a sample in order to show it with grey fill instead of a colored line.")),
      downloadButton("downloadData", label = "Download sample data (zipped)"),
      tags$br(),
      tags$p(paste("code version: ", script.version, sep="")),
      tags$a(href="license.pdf", target="_blank", "usage licence"),
      hr(),
      fileInput('file1', 'Load ladder Data (grey fill)', accept='.csv'),
      fileInput('file2', 'Load Sample Data (blue line)', accept='.csv'),
      fileInput('file3', 'Load Sample2 Data (optional, red line)', accept='.csv'),
      hr()
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("plot1")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$downloadData <- downloadHandler(
    filename <- function() { "RBioanalyzer_sample.zip" },
    content <- function(file) {
      file.copy("www/RBioanalyzer_sample.zip", file)
    },
    contentType = "application/zip"
  )
  
  ladder <- reactive({
    # load ladder data
    inFile1 <- input$file1
    if (is.null(inFile1)) return(NULL)
    ladder <- read_csv(inFile1$datapath, skip = 17)
    ladder <- head(ladder, -1)
    data.ladder <- data.frame("Time"=as.numeric(ladder$Time), "Value"=as.numeric(ladder$Value))
    data.ladder$Value[data.ladder$Value < 0] <- 0
    data.ladder
  })
  
  sample <- reactive({
    # load sample data
    inFile2 <- input$file2
    if (is.null(inFile2)) return(NULL)
    smpl <- read_csv(inFile2$datapath, skip = 17)
    smpl <- head(smpl, -1)
    data.smpl <- data.frame("Time"=as.numeric(smpl$Time), "Value"=as.numeric(smpl$Value))
    data.smpl$Value[data.smpl$Value < 0] <- 0
    data.smpl
  })
  
  sample2 <- reactive({
    # load sample2 data
    inFile3 <- input$file3
    if (is.null(inFile3)) return(NULL)
    smpl2 <- read_csv(inFile3$datapath, skip = 17)
    smpl2 <- head(smpl2, -1)
    data.smpl2 <- data.frame("Time"=as.numeric(smpl2$Time), "Value"=as.numeric(smpl2$Value))
    data.smpl2$Value[data.smpl2$Value < 0] <- 0
    data.smpl2
  })
  
  output$plot1 <- renderPlot({
    if ((is.null(ladder()) || is.null(sample()))) return(NULL)
    
    # plot sample2
    if (is.null(sample2())){
      foo_key <- list(x = .97, y = .92, corner = c(1, 1),
                      text = list(c(input$file1$name, input$file2$name), col = c("grey", "blue")),
                      lines = list(type = c("l", "l"), col = c("grey", "blue"), lwd = 2))
      # plot ladder
      a <- lattice::xyplot(Value ~ Time, ladder(), 
                  panel = function(x,y, ...){
                    panel.xyplot(x,y, type = "l", lwd=0.1, col.line = 'white')
                    panel.polygon(x,y, alpha=0.25, ..., col='grey')
                  })
      
      # plot sample1
      b <- lattice::xyplot(Value ~ Time, sample(), 
                  type = "l", lwd=2, col.line = "blue",
                  key=foo_key,
                  panel = function(x,y, ...){
                    panel.xyplot(x,y, ...)
                  })
      
      b + latticeExtra::as.layer(a)
    } else {
      foo_key2 <- list(x = .97, y = .92, corner = c(1, 1),
                       text = list(c(input$file1$name, input$file2$name, input$file3$name), col = c("grey", "blue", "red")),
                       lines = list(type = c("l", "l", "l"), col = c("grey", "blue", "red"), lwd = 2))
      # plot ladder
      a <- lattice::xyplot(Value ~ Time, ladder(), 
                  panel = function(x,y, ...){
                    panel.xyplot(x,y, type = "l", lwd=0.1, col.line = 'white')
                    panel.polygon(x,y, alpha=0.25, ..., col='grey')
                  })
      
      # plot sample1
      b <- lattice::xyplot(Value ~ Time, sample(), 
                  type = "l", lwd=2, col.line = "blue",
                  key=foo_key2,
                  panel = function(x,y, ...){
                    panel.xyplot(x,y, ...)
                  })
      
      # plot sample2
      c <- lattice::xyplot(Value ~ Time, sample2(), 
                  type = "l", lwd=2, col.line = "red",
                  key = foo_key2,
                  panel = function(x,y, ...){
                    panel.xyplot(x,y, ...)
                  })
      
      b + latticeExtra::as.layer(c) + latticeExtra::as.layer(a)
    }
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

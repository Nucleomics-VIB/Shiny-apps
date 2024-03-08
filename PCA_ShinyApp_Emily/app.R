library(rlist)
#library(devtools)
#devtools::install_github("vqv/ggbiplot")
library(ggbiplot)
library(shinyjs)
library(shinyBS)
library(ggplot2)
library(ggrepel)
library(grDevices)
library(EDASeq)

# name to display on top
app.name <- "PCA plot of RNAseq normalized counts"

# allow larger uploads
options(shiny.maxRequestSize=1000*1024^2)

#######
# ui.r
#######

ui <- fluidPage(
  
  shinyjs::useShinyjs(),
  # App title ----
  #titlePanel("Uploading Files"),
  
  # Application header
  headerPanel("PCA plot of RNAseq normalized counts"),
  windowTitle = "RNASeq PCA plot",
  tags$a(href="https://nucleomicscore.sites.vib.be/en", 
      target="_blank",
      img(src="nucleomics_core_cmyk_pos.jpg", 
          align = "right", 
          width="150", 
          height="58.5", 
          alt="VIB Nucleomics Core")
      ),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      print("You will have to email nucleomics.bioinformatics@vib.be to access the .RData file and metadata .txt file for your project."),
      hr(),
      # download example data
      downloadButton("downloadRData", label = "Download sample RData"),
      downloadButton("downloadMeta", label = "Download sample metadata file"),
      hr(),
      # upload user data
      tipify(fileInput("file1", "Choose RData File",
                multiple = FALSE,
                accept = c(".RData")), 
                "This file will contain the normalized RNA-seq counts per sample and per gene that were processed using the newSeqExpressionSet() function of EDASeq. This data will be processed to create the PCA plots."),
      tipify(fileInput("file2", "Choose metadata .txt File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".txt/.csv")),
                "This file will contain the sample name and group information. If you would like to exclude outliers from your PCA, please delete them from this .txt file and then upload here."),
      textInput("filename","Set filename for .png/.eps file", value = "PCAplot"),
      textInput("MMfile","Set filename for Materials & Methods", value = "M&M"),
      radioButtons("ellipses", "Show circle around groups?",
                   c("Yes" = TRUE,
                     "No" = FALSE), 
                   selected = TRUE),
      radioButtons("labels", "Show labels for sample names?",
                   c("Yes" = TRUE,
                     "No" = FALSE), 
                   selected = TRUE),
      disabled(downloadButton("downloadPCAplot", "Download PCA .png file")),
      disabled(downloadButton("downloadEPS", "Download PCA .eps file")),
      tipify(downloadButton('downloadMM', 'Download M&M'),"Download a text including the names and versions of all packages used in this webtool"),
      hr(),
      print("Created by Emily Sheridan")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      # show plot
      plotOutput("plotPCA", height = "800px", width = "800px"),
      br(),
      # show PCA stats
      tableOutput("pcastats")
      )
    
  )
)

###########
# server.r
###########

server <- function(input, output, session) {
  
  output$downloadRData <- downloadHandler(
    filename <- function() { paste("sample-PpData-Excel", "RData", sep=".") },
    content <- function(file) { file.copy("www/sample-PpData-Excel.RData", file) }
  )
  
  output$downloadMeta <- downloadHandler(
    filename <- function() { paste("sample-metadata", "txt", sep=".") },
    content <- function(file) { file.copy("www/sample-metadata.txt", file) }
  )
  
#  fileOptions <- reactiveValues(currentOptions=c("D.B."))
  
#   observeEvent(input$file1, {
#     fileOptions$currentOptions = list.append(fileOptions$currentOptions, input$file1$datapath)
#   })
#   
#   observeEvent(input$file2, {
#     fileOptions$currentOptions = list.append(fileOptions$currentOptions, input$file2$datapath)
#   })
  
  # import RData object con,taining the counts
  loadData <- reactive({
    req(input$file1)
    # Use a reactiveFileReader to read the file on change, and load the content into a new environment
    load(input$file1$datapath, envir = .GlobalEnv)
    PpData.forExcel
  })
  
  # load metadata file with sample names and group names
  metadata <- reactive({
    req(input$file2)
    sample.info <- read.table(file=input$file2$datapath,sep="\t",quote="",
                              col.names=c("Sample","Group"),
                              header=F,
                              #row.names = 1,
                              as.is=T,stringsAsFactors = T)
    # sample.info$Sample.dot <- gsub("@",".", sample.info$Sample)
    # Group <- factor(sample.info$Group) # used later
    sample.info
  })

  # compute PCA from the loaded data
  computePCA <- reactive({
    #req(input$file2)
    if (is.null(loadData())) return(NULL)
    temp <- t(normCounts(loadData()))
    # subset samples to those included in teh metadata uploaded file
    # this allows uploading an edited list of samples and omit some
    temp <- temp[metadata()$Sample,]
    counts.pca <- prcomp(temp[,apply(temp,2,sd)!=0], scale. = TRUE)
  })
  
  sample.labels <- function()({
    # Split the sample names on '@' and remove last 3 elements, then paste using '_'
    my_delimiter <- '_'
    simpfun <- function(x='', delimiter='_') {
      lst <- strsplit(x, "@")[[1]]
      lst2 <- lst[1:(length(lst) - 3)]
      label <- paste(lst2, collapse = delimiter)
      return(label)
      }
    labels <- lapply(rownames(computePCA()$x), FUN=simpfun, delimiter=my_delimiter)
    labels
    })
  
  # show or not ellipses
  drawEllipses <- function(){
    if (is.null(computePCA())) return(NULL)
    if (length(rownames(computePCA()$x)) > 4){
      drawEllipses <- input$ellipses
    }else{
      drawEllipses <- FALSE
    }
  }
  
  # create a ggplot objet to plot
  preparePlot <- function(){
    if (is.null(computePCA())) return(NULL)
   
    # size of text labels
    labelsize <- 5
    pointsize <- 6
    
    # label transparency, <1 is more tansparent
    alpha <- 1
    
    # other sizes in plot
    My_Theme = theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      legend.text = element_text(size = 12))
    
    # base plot for PCA2 vs PCA1 (choices)
    g <- ggbiplot(computePCA(), 
                  choices = 1:2,
                  obs.scale = 1, 
                  var.scale = 1,
                  var.axes = F,
                  groups = metadata()$Group,
                  ellipse = as.logical(drawEllipses()), 
                  circle = TRUE,
                  ellipse.prob = 0.68,
    )
    g <- g + geom_point(aes(colour=metadata()$Group), size = pointsize, alpha=alpha)
    g <- g + scale_color_discrete(name = '')
    g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
    if (length(levels(metadata()$Group))>7){
      g <- g + guides(col = guide_legend(ncol = 5, byrow = TRUE))
    }
    
    # add labels if TRUE
    if (input$labels == TRUE){
      # add labels using grepel
      nudgex <- 7
      nudgey <- 7
      g <- g + geom_text_repel(aes(label = sample.labels()), 
                               size = labelsize, 
                               max.overlaps=10, 
                               nudge_x=nudgex, 
                               nudge_y=nudgey,
                               segment.colour = NA)
      }
      else{
        # nothing to add (yet!)
      }
    
    # set theme and output plot
    g + My_Theme
  }
  
  # display PCA plot
  output$plotPCA <- renderPlot({
    if (is.null(preparePlot())) return(NULL)
    print(preparePlot())
    shinyjs::enable("downloadPCAplot") # allows the download plot buttons to become active
    shinyjs::enable("downloadEPS")
  })
  
  # display PCA data below the plot
  output$pcastats <- renderTable({
    if (is.null(computePCA())) return(NULL)
    pca.results <- data.frame(summary(computePCA())$importance)
    pca.results[1,] <- round(pca.results[1,], 2)
    pca.results[c(2,3),] <- lapply(pca.results[c(2,3),], function(x) sprintf("%0.1f%%", x*100.0))
    pca.results
    }, 
    rownames = TRUE)
  outputpng <- function(){
    if (is.null(computePCA())) return(NULL)
    png(file = computePCA(),
        units="cm",
        width = 30, height = 30,
        res = 300, bg = "white")
  }
  
  output$downloadPCAplot <- downloadHandler(
    # if (is.null(preparePlot())) return(NULL)
    filename <- function(){paste0(input$filename, ".png")},
    content <- function(file){
      ggsave(file, plot = preparePlot(), device = "png",
             units="cm",
             width = 32, height = 30)
     }
  )
  
  output$downloadEPS <- downloadHandler(
    filename <- function(){paste0(input$filename, ".eps")},
    content <- function(file){
    ggsave(file, plot = preparePlot(), device = "eps")
    }
  )

  output$downloadMM <- downloadHandler(
    filename <- function(){
      paste(input$MMfile, "_session_info.txt", sep="")
    },
    content = function(file) {
      sink(file, append=TRUE)
      cat(paste("Thanks for using our tool:", app.name, "\n", sep=" "))
      cat("\nYou can contact VIB Nucleomics Core at nucleomics@vib.be for any question\n")
      cat(paste("This data was generated on ", format(Sys.time(), "%a %b %d %H:%M:%S %Y"), "\n",sep=" "))
      cat(paste("\nUsing the Principal Component Analysis (PCA), we can project the samples on a two-dimensional graph using the
two first principal components that explain the best the biological variation between those samples.
Each point corresponds to a sample plotted by Principal Component 1 (PC1) and Principal Component 2 (PC2). The ellipses are 68% data ellipses for each
of the groups of samples in the data. The plot can be used to examine the samples for outliers and other
relationships. When normalization successfully removed technical artefacts, the relative distances should be
biologically interpretable.\n"))
      cat("\nThe R packages used in the tools are listed next:\n")
      print(capture.output(sessionInfo()))
      sink()
    }
  )
  
}

######################
# Run the application
######################

shinyApp(ui = ui, server = server)

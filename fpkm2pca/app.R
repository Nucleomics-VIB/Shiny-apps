# fpkm2pca.shinyapp
# A R/shiny tool to create a correlation and PCA plot plots
# from a Nucleomics Core MS-Excel RNASeq count file

library("shiny")
library("shinyBS")
library("shinybusy")
library("shinyjs")
library("shinyalert")
library("openxlsx")
#library("devtools")
library("ggbiplot")
# install_github("vqv/ggbiplot")
#library("ggcorrplot")
# devtools::install_github("kassambara/ggcorrplot")
library("gplots")
library("readr")
library("marray")
library("RColorBrewer")        
library("DT")
#library("pheatmap")

# you may uncomment the next line to allow large input files
options(shiny.maxRequestSize=1000*1024^2)
# the following test checks if we are running on shinnyapps.io to limit file size dynamically
# ref: https://stackoverflow.com/questions/31423144/how-to-know-if-the-app-is-running-at-local-or-on-server-r-shiny/31425801#31425801
#if ( Sys.getenv('SHINY_PORT') == "" ) { options(shiny.maxRequestSize=1000*1024^2) }

app.name <- "fpkm2pca"
script.version <- "1.0"

# Define UI for application that draws a histogram
ui <- fluidPage(

  shinyjs::useShinyjs(),
  #add_busy_gif(
  #  src = "https://jeroen.github.io/images/banana.gif",
  #  height = 70, width = 70
  #),
  #add_busy_spinner(spin = "fading-circle", timeout = 1000),
  add_busy_spinner(spin = "atom", timeout = 1000),
  
  # Application header
  headerPanel("Correlation and PCA plots
              from RNASeq fpkm data"),
  
  # Application title
  titlePanel(
    windowTitle = "RNASeq variance analysis (PCA)",
    tags$a(href="https://corefacilities.vib.be/nc", target="_blank",
           img(src='logo.png', align = "right", width="150", height="58.5", alt="VIB Nucleomics Core"))
  ),
  
  # Sidebar with input
  sidebarLayout(
    # show file import and molecule filters
    sidebarPanel(
      tags$h5(paste(app.name, " version: ", script.version, sep="")),
      downloadButton("downloadCountData", label = "Download test data"),
      downloadButton("downloadSampleGroups", label = "Download test groups"),
      tags$br(),
      tags$a(href="license.pdf", target="_blank", "usage licence"),
      tags$br(),
      tags$a(href="javascript:history.go(0)", tags$i("reset page content"), alt="Reset page"),
      tags$hr(),
      tipify(fileInput('file1', 'Choose RNASeq XLSX count File', accept='.xlsx'),
             "the Data is a MS-Excel file provided by the Nucleomics Core, with worksheet#2 reporting gene expression (FPKM), you may produce a compatible file based on the test data provided here."),
      tipify(fileInput('file2', 'Choose Sample Group File', accept='.txt'),
             "the sample group file is a two-columns text files containing a column of sample names and a column with matching group names"),
      tags$hr(),
      tags$h4("Edit settings and click: ", tags$em("Process & Plot")),
      actionButton(inputId='goButton', "Process & Plot", style='padding:4px; font-weight: bold; font-size:150%'),
      tags$br(),
      tags$hr(),
      textInput('outfile', "name prefix for output files:", value="plot"),
      tipify(radioButtons(inputId = "pcax", 
                          label = "Principal component for X:",
                          choices = 1:4,
                          selected = "1",
                          inline = "TRUE"),
             "choose a principal componenty for the X-axis"),
      tipify(radioButtons(inputId = "pcay", 
                          label = "Principal component for Y:",
                          choices = 1:4,
                          selected = "2",
                          inline = "TRUE"),
             "choose a principal componenty for the Y-axis"),
      tipify(radioButtons(inputId = "ellipse", 
                          label = "Draw ellipses:",
                          choices = c("Yes" = TRUE, "No" = FALSE),
                          selected = TRUE,
                          inline = "TRUE"),
             "whether to draw an ellipse around each group"),
      tipify(radioButtons(inputId = "showlabel", 
                          label = "Show labels:",
                          choices = c("Yes" = TRUE, "No" = FALSE),
                          selected = TRUE,
                          inline = "TRUE"),
             "whether to show sample labels next to the points"),
      tipify(sliderInput("size", "label size:",
                                 min = 1, max = 10,
                                 value = 3, step = 1,
                                 animate = TRUE),
             "choose a font / label size"),
      tipify(downloadButton('downloadMM', 'Download M&M'),"Download a text including the names and versions of all packages used in this webtool")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      
      tabsetPanel(type = "tabs",
                  tabPanel("Info", 
                           h3("RNASeq Data"),
                           h4("The results of your RNASeq sequencing typically include two main Excel files:"),
                           tags$div(tags$ul(
                             tags$li("a count file returning raw and normalise gene counts (used here)"),
                             tags$li("a statistical analysis file returning pairwise differential expression values and statistics obtained after comparing sample groups."),
                           ), style = "font-size: 15px"),
                           p("If not provided by us, you will also need a tab-separated 'Sample Group file' with two columns"),
                           tags$div(tags$ul(
                             tags$li("one containing the sample names (available as column names in the count file)"),
                             tags$li("the second containing the groups to which these samples belong."),
                           ), style = "font-size: 15px"),                           
                           h3("Variance analysis"),
                           p("One way to estimate the quality of your experiment is to confirm that your sample groups are mutually distinct at the level of gene expression."),
                           h4("We perform here such analysis using two methods:"),
                           p("- plotting the Spearman-correlation values of all-versus-all sample comparisons based on the normalized counts (Correlation Plot)"),
                           p("- plotting pairs of principal components (PCs) after computing a principal component analysis (PCA) of the normalized counts (PCA plot)"),
                           h3("Interpretation"),
                           p("If your experiment went well and your conditions induce differential expression of genes, you should see your samples grouped together by group in the heatmap and the groups clearly separated in the PCA plots."),
                           p("If a sample is shown in a different group, this may be due to experimental errors (eg. sample swap) when the group are furthermore nicely distinct."),
                           p("Note: when drawn, the ellipse shows the 68% probability limit for each group."),
                           p("By contrast, overlapping PCA clouds may indicate that the conditions applied to the samples did not induce enough differential expression between groups to separate them based on variance."),
                           br(),
                           p("Please contact us at nucleomics.bioinformatics@vib.be for more information about how to interpret these plots or if you encounter problems using this tool."),
                           hr(),
                           textOutput('full.data.cnt'),
                           textOutput('filt.data.cnt'),
                  ),
                  tabPanel("Correlation Plot", 
                           plotOutput("corrplot", width = "80%", height = "600px"),
                           disabled(downloadButton('downloadCorrPlot', 'Download Plot'))
                  ),
                  tabPanel("PCA Plot", 
                           plotOutput("pcaplot", width = "80%", height = "600px"),
                           disabled(downloadButton('downloadPCAPlot', 'Download Plot'))
                  )
      )
    )
  )
  # end UI block
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$showellipse.value <- renderText({ input$ellipse })  
  output$showlabel.value <- renderText({ input$showlabel })
  output$pcax.value <- renderText({ input$pcax })
  output$pcay.value <- renderText({ input$pcay })
  
  # download test count data file
  output$downloadCountData <- downloadHandler(
    filename <- function() { paste("expXXXX-RNAseqCounts", "xlsx", sep=".") },
    content <- function(file) { file.copy("www/expXXXX-RNAseqCounts.xlsx", file) }
  )
  
  # download test sample group file 
  output$downloadSampleGroups <- downloadHandler(
    filename <- function() { paste("XXXX-sample-groups", "txt", sep=".") },
    content <- function(file) { file.copy("www/XXXX-sample-groups.txt", file) }
  )
  
  # import user provided count data
  fpkm.data <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    
    # load data from excel file
    dat <- read.xlsx(inFile$datapath, sheet=2)
    
    # keep only data columns (remove last columns starting from "Chromosome")
    chromosome.col <- which(colnames(dat)==as.vector("Chromosome"))
    fpkm.data <- dat[,c(3:(chromosome.col-1))]
    
    # simplify sample names
    colnames(fpkm.data) <- sapply(colnames(fpkm.data), function(strings){
      ind = unlist(gregexpr(pattern = "@", text = strings))
      if (length(ind) < 3){NA}
      else{substr(strings, 1, ind[length(ind) - 2] - 1)}
    })
    
    # kick useless part of names for samples
    colnames(fpkm.data) <- sub("@", "_", colnames(fpkm.data))

    # return data as 'fpkm.data()'
    fpkm.data
  })
  
  # show full data row count
  output$full.data.cnt <- reactive({
    #if (is.null(fpkm.data())) return("Waiting for data!")
    paste("Rows in the Full data: ", nrow(fpkm.data()))
  })
  
  # import user provided sample group file
  sample.groups <- reactive({
    inFile <- input$file2
    if (is.null(inFile)) return(NULL)
    
    sample.groups <- read_delim(inFile$datapath, 
                                delim = "\t", 
                                escape_double = FALSE, 
                                col_names = FALSE, 
                                trim_ws = TRUE,
                                show_col_types = FALSE)
    
    colnames(sample.groups) <- c("labels", "group")
    
    # simplify sample labels
    sample.groups$labels <- sapply(sample.groups$labels, function(strings){
      ind = unlist(gregexpr(pattern = "@", text = strings))
      if (length(ind) < 3){NA}
      else{substr(strings, 1, ind[length(ind) - 2] - 1)}
    })
    
    # kick useless part of names for samples
    sample.groups$labels <- sub("@", "_", sample.groups$labels)

    # return sample.groups()'
    sample.groups
  })
  
  # process data  
  filtered.data <- eventReactive({input$goButton}, {
    # do nothing in absence of data
    if (is.null(fpkm.data())) return(NULL)
    if (is.null(sample.groups())) return(NULL)
    
    # select only rows with variance
    counts <- as.data.frame(fpkm.data())
    filtered.data <- counts[apply(counts, 1, var) != 0,]
    
    # return as filtered.data()
    filtered.data
  })
  
  # show filtered data row counts
  output$filt.data.cnt <- reactive({
    if (is.null(filtered.data())) return("Waiting for data!")
    paste("Rows in the filtered data: ", nrow(filtered.data()))
  })
  
  # prepare Corr plot
  corrplot.data <- function(){
    if (is.null(filtered.data())) return(NULL)
    
    counts.cor <- cor(filtered.data(), use="pairwise.complete.obs", method="spearman")
    # create palette with the marray package
    my.pal <- marray::maPalette(low="green", high="red", mid="yellow", k=69)
    par(mar=c(7,4,4,2)+0.1)
    corrplot.data <- heatmap.2(counts.cor,
                               col=my.pal,
                               cexRow=1,
                               cexCol=1,
                               trace="none",
                               scale="none",
                               margins=c(8,8))
    corrplot.data
  }
  
  # show Corr plot on page
  output$corrplot <- renderPlot({
    if (is.null(corrplot.data())) return(NULL)
    par(mar=c(7,4,4,2)+0.1)
    print(corrplot.data())
    shinyjs::enable("downloadCorrPlot")
  })
  
  # download Corr plot
  output$downloadCorrPlot2 <- downloadHandler(
    if (is.null(corrplot.data())) return(NULL),
    filename = function() { 
      paste0(input$outfile, "_Corr.png", sep='') 
      },
    content = function(file) {
      png(file)
      print(corrplot.data())
      null <- dev.off()
      }
  )

  # prepare Corr plot
  pcaplot.data <- function(){
    if (is.null(filtered.data())) return(NULL)
    if (as.integer(input$pcax) == as.integer(input$pcay)) {
      shinyalert("Oops!", "PCX should be different of PCY", type = "error")
      return(NULL)
      }
    counts.pca <- prcomp(t(filtered.data()), scale. = TRUE)
    pca.table <- as.data.frame(counts.pca$x[,1:4])
    pca.table$labels <- rownames(pca.table)
    # ggbiplot can only draw ellipses if there are at least 5 samples
    drawEllipses <- ifelse(length(rownames(counts.pca$x)) > 4, input$ellipse, FALSE)
    
    # create plot
    # show labels YES or NO from user choice
    if ( input$showlabel == "TRUE" ) {
      g <- ggbiplot(counts.pca, 
                  choices = c(as.integer(input$pcax), as.integer(input$pcay)), 
                  obs.scale = 1, 
                  var.scale = 1, 
                  var.axes = F,
                  groups = sample.groups()$group,
                  ellipse = as.logical(drawEllipses), 
                  ellipse.prob = 0.68,
                  circle = TRUE,
                  labels = sample.groups()$labels,
                  labels.size = as.integer(input$size)
                  )
    } else {
      g <- ggbiplot(counts.pca, 
                    choices = c(as.integer(input$pcax), as.integer(input$pcay)), 
                    obs.scale = 1, 
                    var.scale = 1, 
                    var.axes = F,
                    groups = sample.groups()$group,
                    ellipse = as.logical(drawEllipses), 
                    ellipse.prob = 0.68,
                    circle = TRUE
      ) + geom_point(aes(colour=sample.groups()$group), size = as.integer(input$size))
    }
    # tune plot format
    g <- g + scale_color_discrete(name = '')
    g <- g + theme(legend.direction = 'horizontal',
                   legend.position = 'top',
                   legend.text=element_text(size=as.integer(3*input$size)),
                   axis.title.x=element_text(size=3*input$size),
                   axis.title.y=element_text(size=3*input$size))
    if (length(unique(pca.table$group))>7){
      g <- g + guides(col = guide_legend(ncol = 5, byrow = TRUE))
    }
    g
  }
  
  # show PCA plot on page
  output$pcaplot <- renderPlot({
    if (is.null(pcaplot.data())) return(NULL)
    print(pcaplot.data())
    shinyjs::enable("downloadPCAPlot")
  })
  
  # download PCA plot
  output$downloadPCAPlot2 <- downloadHandler(
      if (is.null(pcaplot.data())) return(NULL),
      filename = function() { 
        paste0(input$outfile, "_PCA.png", sep='') 
      },
      content = function(file) {
        png(file)
        print(pcaplot.data())
        null <- dev.off()
      }
    )
  
  # download M&M text
  output$downloadMM <- downloadHandler(
    filename = function() {
      paste(input$outfile, "_session_info.txt", sep="")
    },
    content = function(file) {
      sink(file, append=TRUE)
      cat(paste("Thanks for using our tool", app.name, script.version, "\n", sep=" "))
      cat ("This tool generates a Heatmap after correlating samples based on their FPKM values and a PCA plot based on the same data.")
      cat("\nYou can contact The Nucleomics Core at nucleomics.bioinformatics@vib.be with your question\n")
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

# counts2pca.shinyapp
# A R/shiny tool to create a correlation and PCA plot plots
# from a Nucleomics Core MS-Excel RNASeq count file
# SP@NC, 2023-01-06

library("shiny")
library("shinyBS")
library("shinybusy")
library("shinyjs")
library("shinyalert")
library("openxlsx")
#library("devtools")
# devtools::install_github("vqv/ggbiplot")
library("ggbiplot")
# devtools::install_github("kassambara/ggcorrplot")
#library("ggcorrplot")
library("gplots")
library("readr")
library("marray")
library("RColorBrewer")
library("DT")
#library("pheatmap")
library("factoextra")

# you may uncomment the next line to allow large input files
options(shiny.maxRequestSize=1000*1024^2)
# the following test checks if we are running on shinnyapps.io to limit file size dynamically
# ref: https://stackoverflow.com/questions/31423144/how-to-know-if-the-app-is-running-at-local-or-on-server-r-shiny/31425801#31425801
#if ( Sys.getenv('SHINY_PORT') == "" ) { options(shiny.maxRequestSize=1000*1024^2) }

app.name <- "counts2pca"
script.version <- "1.3b"

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
              from RNASeq raw count data"),
  
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
             "the Data is a MS-Excel file provided by the Nucleomics Core, with worksheet#1 reporting gene expression counts, you may produce a compatible file based on the test data provided here."),
      tipify(fileInput('file2', 'Choose Sample Group File', accept='.txt'),
             "the Sample Group file is a tab-separated text files with a header row (more info on the right)"),
      tipify(numericInput("mincnt", 
                          "min SUM(counts):", 
                          1,
                          min = 1, 
                          max = 1000),
             "require a minimum for the sum of all sample counts to keep a gene row"),
      tags$h4("Edit settings above and click: ", tags$em("Process & Plot")),
      actionButton(inputId='goButton', "Process & Plot", style='padding:4px; font-weight: bold; font-size:150%'),
      tags$hr(),
      tags$h4("Plot file name"),
      textInput('outfile', "name prefix for output files:", value="plot"),
      tags$h4("PCA settings (auto-refresh)"),
      tipify(radioButtons(inputId = "pcax", 
                          label = "Principal component for X:",
                          choices = 1:6,
                          selected = "1",
                          inline = "TRUE"),
             "choose a principal componenty for the X-axis"),
      tipify(radioButtons(inputId = "pcay", 
                          label = "Principal component for Y:",
                          choices = 1:6,
                          selected = "2",
                          inline = "TRUE"),
             "choose a principal componenty for the Y-axis"),
      tipify(selectInput(inputId = "groups", 
                         label = "Choose the grouping variable:",
                         choices = c("group"),
                         selected="group"),
             "choose a sample grouping column"),
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
      tipify(downloadButton('downloadMM', 'Download M&M'),
              "Download a text including the names and versions of all packages used in this webtool")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      
      tabsetPanel(type = "tabs",
                  tabPanel("Info", 
                           h3("RNASeq Data"),
                           h4("the results of your RNASeq sequencing typically include two main Excel files:"),
                           tags$div(tags$ul(
                             tags$li("a count file returning raw (used here), and normalise gene counts"),
                             tags$li("a statistical analysis file returning pairwise differential expression values and statistics obtained after comparing sample groups (not used here)."),
                           ), style = "font-size: 15px"),
                           p("If not provided by us, you will also need a tab-separated 'Sample Group' file with a header row and two or more columns"),
                           tags$div(tags$ul(
                             tags$li("the first column should be named 'labels' and contain the sample identifiers found as column names in the Excel file"),
                             tags$li("the second column should be named 'groups' and link each sample to a user-defined group for ellipse plotting in the PCA"),
                             tags$li("optionally: more metadata columns can be added to allow alternative grouping of the PCA results with ellipses"),
                             tags$li("Note: if you want to omit samples from the analysis, simply omit the corresponding rows in the group file."),
                           ), style = "font-size: 15px"),                           
                           h3("Variance analysis"),
                           p("One way to estimate the quality of your experiment is to confirm that your sample groups are mutually distinct at the level of gene expression."),
                           h4("We perform here such analysis using two methods:"),
                           p("- plotting the Spearman-correlation values of all-versus-all sample comparisons based on the raw counts (Correlation Plot)"),
                           p("- plotting pairs of principal components (PCs) after computing a principal component analysis (PCA) of the filtered* raw counts (PCA plot)"),
                           p("(*) filtering removes rows with no variance and rows where the sum of all samples counts is less than the user provided minimum"),
                           h3("Interpretation"),
                           p("If your experiment went well and your conditions induce differential expression of genes, you should see your samples grouped together by group in the correlation heatmap and the groups clearly separated in the PCA plots."),
                           p("If a sample is shown in a different group, this may be due to experimental errors (eg. sample swap) when the group are furthermore nicely distinct."),
                           p("Note: when drawn, the ellipse shows the 68% probability limit for each group (default value in the R  package)"),
                           p("By contrast, overlapping PCA clouds may indicate that the conditions applied to the samples did not induce enough differential expression between groups to separate them based on variance."),
                           br(),
                           p("Please contact us at 'nucleomics.bioinformatics@vib.be' for more information about how to interpret these plots or if you encounter problems using this tool."),
                           hr(),
                           textOutput('full.data.cnt'),
                           textOutput('filt.data.cnt'),
                           textOutput('grouping.variable'),
                           textOutput('grouping.values')
                  ),
                  tabPanel("Correlation Plot", 
                           plotOutput("corrplot", width = "80%", height = "600px"),
                           disabled(downloadButton('downloadCorrPlot', 'Download Plot'))
                  ),
                  tabPanel("PCA Plot",
                           br(),
                           p("the first plot returns the % of variance explained by the first max=6 dimensions of the PCA"),
                           plotOutput("pcavar", width = "50%", height = "250px"),
                           disabled(downloadButton('downloadPCAVar', 'Download Plot')),
                           p("the next plot returns a two-dimension plot based on the user choices"),
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
  count.data <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    
    # load data from excel file (sheet#1 contains raw counts)
    # dat <- read.xlsx("expXXXX-RNAseqCounts.xlsx", sheet=1)
    dat <- read.xlsx(inFile$datapath, sheet=1)
    
    # keep only data columns (remove last columns starting from "Chromosome")
    chromosome.col <- which(colnames(dat)==as.vector("Chromosome"))
    count.data <- dat[,c(3:(chromosome.col-1))]
    
    # simplify sample names
    # colnames(count.data) <- sapply(colnames(count.data), function(strings){
    #   ind = unlist(gregexpr(pattern = "@", text = strings))
    #   if (length(ind) < 3){NA}
    #   else{substr(strings, 1, ind[length(ind) - 2] - 1)}
    # })
    
    # keep only first string in long @-delimited names
    colnames(count.data) <- sapply(colnames(count.data), function(strings){
      ind = unlist(gregexpr(pattern = "@", text = strings))
      if (length(ind) < 3){NA}
      else{substr(strings, 1, ind[1]-1)}
    })
    
    # remove leftover '@' symbols
    colnames(count.data) <- sub("@", "_", colnames(count.data))

    # return data as 'count.data()'
    count.data
  })
  
  # show full data row count
  output$full.data.cnt <- reactive({
    #if (is.null(count.data())) return("Waiting for data!")
    paste("Rows in the Full data: ", nrow(count.data()))
  })
  
  # import user provided sample group file
  sample.groups <- reactive({
    inFile <- input$file2
    if (is.null(inFile)) return(NULL)

    # sample.groups <- read_delim("XXXX-sample-groups.txt",
    #                             delim = "\t",
    #                             escape_double = FALSE,
    #                             col_names = TRUE,
    #                             trim_ws = TRUE,
    #                             show_col_types = FALSE)
    
    sample.groups <- read_delim(inFile$datapath, 
                                delim = "\t", 
                                escape_double = FALSE, 
                                col_names = TRUE, 
                                trim_ws = TRUE,
                                show_col_types = FALSE)

    # check if table row contains the required 2 columns
    if (! all(c("labels", "group") %in% colnames(sample.groups)) ) {
      shinyalert("Oops!", "the group file should have a header row with at least two columns named 'labels' and 'group'", type = "error")
      return(NULL)
    }
    
    # simplify sample labels
    # sample.groups$labels <- sapply(sample.groups$labels, function(strings){
    #   ind = unlist(gregexpr(pattern = "@", text = strings))
    #   if (length(ind) < 3){NA}
    #   else{substr(strings, 1, ind[length(ind) - 2] - 1)}
    # })
    
    # keep only first string in long @-delimited names
    sample.groups$labels <- sapply(sample.groups$labels, function(strings){
         ind = unlist(gregexpr(pattern = "@", text = strings))
         substr(strings, 1, ind[1]-1)
       })
    
    # kick useless part of names for samples
    sample.groups$labels <- sub("@", "_", sample.groups$labels)

    # update the dropdown in the UI based on the extra column names
    items <- as.character(colnames(sample.groups)[-1])
    updateSelectInput(session = getDefaultReactiveDomain(), 
                      inputId = "groups", 
                      choices = items, 
                      selected = "group")
    
    # return sample.groups()
    sample.groups
  })
  
  # filter data  
  filtered.data <- eventReactive({input$goButton}, {
    # do nothing in absence of data
    if (is.null(count.data())) return(NULL)
    if (is.null(sample.groups())) return(NULL)

    # create local data.frame for filtering
    counts <- as.data.frame(count.data())
    
    # keep only samples listed in the group file
    keep.cols <- sample.groups()$labels
    filtered.data <- counts[,keep.cols]
    
    # select only rows with variance
    filtered.data <- filtered.data[apply(filtered.data, 1, var) != 0,]
    
    # keep only rows with sum >= input$mincnt
    keep.rows <- rowSums(filtered.data) >= input$mincnt
    filtered.data <- filtered.data[keep.rows,]
    
    # return as filtered.data()
    filtered.data
  })
  
  # show filtered data row counts
  output$filt.data.cnt <- reactive({
    if (is.null(filtered.data())) return("Waiting for data!")
    paste("Rows in the filtered data: ", nrow(filtered.data()))
  })
  
  # show selected grouping variable
  output$grouping.variable <- reactive({
    if (is.null(input$groups)) return("Waiting for data!")
    paste("Current grouping column: ", input$groups)
  })
  
  # show selected grouping variable
  output$grouping.values <- reactive({
    if (is.null(input$groups) || is.null(sample.groups())) return("Waiting for data!")
    paste("Current grouping values: ", sample.groups()[input$groups])
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
    # par(mar=c(7,4,4,2)+0.1)
    pdf(file = NULL)
    print(corrplot.data())
    shinyjs::enable("downloadCorrPlot")
  })
  
  # download Corr plot
  output$downloadCorrPlot <- downloadHandler(
    if (is.null(corrplot.data())) return(NULL),
    filename = function() { paste(input$outfile, "_Corr.png", sep='') },
    content = function(file) {
      par(mar=c(7,4,4,2)+0.1)
      png(filename = file)
      print(corrplot.data())
      null <- dev.off()
    }
  )

  # prepare PCA plot
  pcaplot.data <- function(){
    if (is.null(filtered.data())) return(NULL)
    if (as.integer(input$pcax) == as.integer(input$pcay)) {
      shinyalert("Oops!", "PCX should be different of PCY", type = "error")
      return(NULL)
    }
    
    # compute PCA from filtered.data
    pca <- prcomp(t(filtered.data()), scale. = TRUE)
    
    # find last dimension in PCA
    dimlim <- ncol(pca$x)
    dim.max <- min(c(dimlim,6))

    # read current values
    pca.x <- as.integer(input$pcax)
    pca.y <- as.integer(input$pcay)
    
    # adapt controls when PCA has less than 6 dimensions
    updateRadioButtons(session = getDefaultReactiveDomain(),
                         inputId = "pcax",
                         choices = 1:dim.max,
                         selected = pca.x,
                         inline = "TRUE"
                       )
    updateRadioButtons(session = getDefaultReactiveDomain(),
                         inputId = "pcay",
                         choices = 1:dim.max,
                         selected = pca.y,
                         inline = "TRUE"
                       )

    # store in table for display
    pca.table <- as.data.frame(pca$x[,1:dim.max])
    pca.table$labels <- rownames(pca.table)
    pca.table$group <- unlist(sample.groups()[input$groups])
    
    # define accompanying class for ellipse grouping
    pca.class <- pca.table$group

    # ggbiplot can only draw ellipses if there are at least 5 samples
    drawEllipses <- ifelse(length(rownames(pca$x)) > 4, input$ellipse, FALSE)

    # create plot
    # show labels YES or NO from user choice
    if ( input$showlabel == "TRUE" ) {
      # show label to TRUE
      
      g <- ggbiplot(pca, 
                    choices = c(pca.x, pca.y), 
                    obs.scale = 1, 
                    var.scale = 1, 
                    var.axes = F,
                    groups = pca.class,
                    ellipse = as.logical(drawEllipses), 
                    ellipse.prob = 0.68,
                    circle = TRUE,
                    labels = sample.groups()$labels,
                    labels.size = as.integer(input$size)
      )
    } else {
      # show label to FALSE
      g <- ggbiplot(pca, 
                    choices = c(pca.x, pca.y), 
                    obs.scale = 1, 
                    var.scale = 1, 
                    var.axes = F,
                    groups = pca.class,
                    ellipse = as.logical(drawEllipses), 
                    ellipse.prob = 0.68,
                    circle = TRUE
      ) + geom_point(aes(colour=pca.class), size = as.integer(input$size))
    }
    # tune plot format
    g <- g + scale_color_discrete(name = '')
    g <- g + theme(legend.direction = 'horizontal',
                   legend.position = 'top',
                   legend.text=element_text(size=as.integer(3*input$size)),
                   axis.title.x=element_text(size=3*input$size),
                   axis.title.y=element_text(size=3*input$size))
    if (length(unique(pca.class))>7){
      g <- g + guides(col = guide_legend(ncol = 5, byrow = TRUE))
    }

    # print variance explained to a second object
    g2 <- fviz_eig(pca, 
                   ncp=6, 
                   main="% variance explained in the first 6 dimensions"
                   )

    # return 2 plots in a list
    return(
      list(plot1=g, plot2=g2)
    )
  }
  
  # show PCA plot on page
  output$pcaplot <- renderPlot({
    if (is.null(pcaplot.data())) return(NULL)
    .tmp <- pcaplot.data()
    print(.tmp$plot1)
    shinyjs::enable("downloadPCAPlot")
  })

  # show PCA var on page
  output$pcavar <- renderPlot({
    if (is.null(pcaplot.data())) return(NULL)
    .tmp <- pcaplot.data()
    print(.tmp$plot2)
    shinyjs::enable("downloadPCAVar")
  })
    
  # download PCA plot
  output$downloadPCAPlot <- downloadHandler(
    if (is.null(pcaplot.data())) return(NULL),
    filename = function() { paste(input$outfile, "_PCA.png", sep='') },
    content = function(file) {
      .tmp <- pcaplot.data()
      ggsave(file, 
             plot = .tmp$plot1, 
             device = "png", 
             width=16, 
             height=16, 
             unit="cm")
    }
  )
  
  # download PCA var
  output$downloadPCAVar <- downloadHandler(
    if (is.null(pcaplot.data())) return(NULL),
    filename = function() { paste(input$outfile, "_PCAvariance.png", sep='') },
    content = function(file) {
      .tmp <- pcaplot.data()
      ggsave(file, 
             plot = .tmp$plot2, 
             device = "png", 
             width=16, 
             height=8, 
             unit="cm")
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
      cat ("This tool generates a PCA plot based on the variance detected in row gene counts from multiple RNASeq samples.")
      cat("\nYou can contact the Nucleomics Core at nucleomics.bioinformatics@vib.be with your question\n")
      cat(paste("This data was generated on ", format(Sys.time(), "%a %b %d %H:%M:%S %Y"), "\n",sep=" "))
      cat("\nthe R packages used in the tools are listed next:\n")
      print(capture.output(sessionInfo()))
      sink()
    }
  )
  
  # end server block
}

# Run the application
shinyApp(ui = ui, server = server)

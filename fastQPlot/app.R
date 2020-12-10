# FastQPlot.shinyapp
# A R/shiny tool to filter fastQ data and make plot
# VIB Nucleomics Core (author: Stephane Plaisance)
# 2020-12-08 v1.0

library("shiny")
library("shinyBS")
library("shinyjs")
library("shinymanager")
library("markdown")
library("ShortRead")
library("DT")
library("stringr")
library("ggplot2")
library("grid")
library("gridExtra")
library("ggpubr")

# App defaults
app.name <- "FastQPlot"
version <- "Version: 1.0 - 2020-12-08"
versnum <- "v1.0"

# load path to Uploads
source("config.R")

# custom functions
textInput3 <- function (inputId, label, value = "", ...) {
  div(
    style = "display:inline-block",
    tags$label(label, `for` = inputId),
    tags$input(
      id = inputId,
      type = "text",
      value = value,
      ...
    )
  )
}

Nvalue <- function(lim, x, na.rm = TRUE) {
  # handle NA values
  if (isTRUE(na.rm)) {
    x <- x[!is.na(x)]
  }
  cutval <- 100 / lim
  # compute LXX and NXX
  sorted <- sort(x, decreasing = TRUE)
  SXX <- sum(x) / cutval
  csum <- cumsum(sorted)
  GTLXX <- as.vector(csum >= SXX)
  LXX = min(which(GTLXX == TRUE))
  NXX <- round(sorted[LXX], 1)
  # eg: get NXX with lst['NXX']
  NXX
}

# you may un-comment the next line to allow 10MB input files
options(shiny.maxRequestSize = 100 * 1024 ^ 2)

ui <- navbarPage(
  paste0(app.name, "! ", versnum),
  id = "tabSwitch",
  tabPanel("Data Upload",
           sidebarLayout(
             sidebarPanel(
               tags$h4(tags$i(
                 "Analyze long read FastQ data (CCS or denovo assembly)"
               )),
               tags$p(
                 "The reads will be subjected to filtering accoring to your chosen limits and the resulting data used to create plots and reveal trends in your data."
               ),
               tags$p(
                 "The filtered reads can be downloaded to your computer (FastQ format) and fed to downstream tools of your choice."
               ),
               tags$p(
                 "Note: The size of the FastQ file should be less than 100Mb. For larger files, we recommend to install this shiny application (and its dependencies) on your own server from our GitHUB copy, see 'Info/About'."
               ),
               tipify(
                 fileInput('FASTQ', 'Upload a FastQ.gz File', accept = '.gz'),
                 "the Data is a bgzipped FastQ file (for the sake of upload speed)"
               ),
               actionButton("importFASTQ", "import and analyze", icon =
                              icon('upload'))
             ),
             mainPanel()
           )),
  tabPanel("Analyze FastQ",
           sidebarLayout(
             sidebarPanel(
               img(src = "fastq.png", style = 'padding-right:10px;height:123px;width:238px'),
               br(),
               tags$h4(tags$i("FastQ filtering")),
               tags$p("Choose limits using the controls below."),
               # read length range
               tags$h5("Read length"),
               numericInput(
                 inputId = "minlen",
                 label = "min",
                 value = 0
               ),
               numericInput(
                 inputId = "maxlen",
                 label = "max",
                 value = 10000
               ),
               tags$h5("Average Base Quality"),
               # minimum
               numericInput("minqual",
                            label = "min",
                            value = 0),
               # maximum
               numericInput("maxqual",
                            label = "max",
                            value = 40),
               tags$h5("GC%"),
               # minimum
               numericInput("minGC",
                            label = "min",
                            value = 0),
               # maximum
               numericInput("maxGC",
                            label = "max",
                            value = 100),
               br(),
               actionButton(
                 "filterFASTQ",
                 "create filtered FastQ",
                 icon =
                   icon('filter', lib = "glyphicon")
               ),
               hr(),
               downloadButton(
                 "downloadFASTQ",
                 "download filtered FastQ",
                 icon =
                   icon('save', lib = "glyphicon")
               ),
               width = 2
             ),
             mainPanel(
               selectInput(
                 inputId = "format",
                 label = "format:",
                 choices = c("PDF" = "pdf", "PNG" = "png"),
                 selected = "png"
               ),
               downloadButton("downloadPlot", "save plot", icon =
                                icon('save', lib = "glyphicon")),
               hr(),
               plotOutput('plots', width = "100%", height = "800px")
             )
           ),),
  navbarMenu(
    "Info",
    tabPanel("About",
             fluidRow(
               column(
                 3,
                 br(),
                 br(),
                 img(src = 'logo.png', width = 200),
                 br(),
                 tags$small(a(href = "https://www.nucleomics.be", "VIB Nucleomics Core"))
               ),
               column(6,
                      includeMarkdown("www/about.md"))
             )),
    tabPanel("License",
             fluidRow(column(
               6,
               includeMarkdown("www/license.md")
             ))),
    tabPanel("Version",
             helpText(version))
  ),
  # initialize Shinyjs: used for click()
  useShinyjs()
  
  # ui ends here
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  # clean old data
  file.remove(file.path(uploads, "data.fq.gz"))
  file.remove(file.path(uploads, "filtered_data.fq.gz"))
  file.remove(file.path(uploads, "data.txt"))
  file.remove(file.path(uploads, "filtered_data.txt"))
  fqdata <- ""
  
  # get user to select a FASTQ file to process
  FASTQFile <- reactive({
    file <- input$FASTQ
    req(file)
    fpath <- file$datapath
    fext <- tools::file_ext(fpath)
    # report to user
    # showModal(modalDialog(
    #     title = "Chosen file",
    #     paste0(fpath, " - " , fext),
    #     easyClose = TRUE,
    #     footer = NULL
    # ))
    validate(need(fext %in% "gz", "Upload must be a .fq.gz or .fastq.gz file"))
    res <-
      list(
        upath = file$datapath,
        fname = basename(file$name),
        fext = fext
      )
    return(res)
  })
  
  # upload FASTQ file
  observeEvent(input$importFASTQ, {
    isolate(FASTQFile()$upath)
    upload_result <- list()
    # copy file to server
    err0 <- file.copy(FASTQFile()$upath, uploads)
    err1 <- file.remove(file.path(uploads, "data.fq.gz"))
    err2 <- file.remove(file.path(uploads, "data.txt"))
    err3 <-
      file.rename(file.path(uploads, paste0("0.", FASTQFile()$fext)),
                  file.path(uploads, "data.fq.gz"))
    if (!err3) {
      upload_result <- c(upload_result, "FAILED! Error during upload")
    }
    # all went well
    if (length(upload_result) == 0) {
      upload_result <-
        c(upload_result, "The file was uploaded succesfully")
    }
    
    # report to user
    showModal(modalDialog(
      title = "Upload status",
      paste(unlist(upload_result), collapse = ", "),
      easyClose = TRUE,
      footer = NULL
    ))
    
    # load fastq data into R object
    fastq.file <- file.path("Uploads", "data.fq.gz")
    f <- FastqFile(fastq.file)
    fqdata <- readFastq(f)
    close(f)
    
    # create data.frame with read metrics
    data <- cbind(
      rname = data.frame(id(fqdata)),
      len = width(fqdata),
      topqual = alphabetScore(fqdata)
    )
    data$avgQual <- round(data$topqual / data$len, 2)
    data <- data[, -3]
    colnames(data) <- c("name", "length", "avgQual")
    
    seq <- data.frame(sread(fqdata))
    seq$G <- str_count(seq$sread.fqdata., "G")
    seq$C <- str_count(seq$sread.fqdata., "C")
    seq$A <- str_count(seq$sread.fqdata., "A")
    seq$T <- str_count(seq$sread.fqdata., "T")
    seq$GC <- round(100 * (seq$G + seq$C) / (seq$A + seq$T + seq$G + seq$C), 2)
    
    data <- cbind(data, seq$GC)
    colnames(data) <- c("name", "length", "avgQual", "GC")
    
    # update NumericInputs from data
    updateNumericInput(session, "minlen", value = min(data$length))
    updateNumericInput(session, "maxlen", value = max(data$length))
    updateNumericInput(session, "minqual", value = min(data$avgQual))
    updateNumericInput(session, "maxqual", value = max(data$avgQual))
    updateNumericInput(session, "minGC", value = min(data$GC))
    updateNumericInput(session, "maxGC", value = max(data$GC))
    
    write.table(
      data,
      file = "Uploads/data.txt",
      row.names = FALSE,
      col.names = TRUE,
      sep = ','
    )
    
    # report to user
    showModal(modalDialog(
      title = "Analysis status",
      paste("analysis done for", nrow(data), "reads", collapse = " "),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  # get full data
  fastqdata <- reactive({
    if (is.null(input$importFASTQ))
      return(NULL)
    
    req("Uploads/data.txt")
    df <- read.csv("Uploads/data.txt",
                   header = TRUE,
                   sep = ",")
    return(df)
  })
  
  # create filtered subset
  #filtered.data <- eventReactive({
  #  input$filter
  #}, {
  filtered.data <- reactive({
    # do nothing in absence of data
    if (is.null(fastqdata()))
      return(NULL)
    fdata <- subset(
      fastqdata(),
      length >= as.numeric(input$minlen) &
        length <= as.numeric(input$maxlen) &
        avgQual >= as.numeric(input$minqual) &
        avgQual <= as.numeric(input$maxqual) &
        GC >= as.numeric(input$minGC) &
        GC <= as.numeric(input$maxGC)
    )
    
    # write to disk
    write.table(
      fdata,
      file = "Uploads/filtered_data.txt",
      row.names = FALSE,
      col.names = TRUE,
      sep = ','
    )
    
    # return filtered table
    return(fdata)
  })
  
  filtered.summary <- reactive({
    # do nothing in absence of data
    if (is.null(filtered.data()))
      return(NULL)
    
    data <- filtered.data()
    
    # get N50's
    n50length <- Nvalue(50, data$length)
    n50qual <- Nvalue(50, data$avgQual)
    n50GC <- Nvalue(50, data$GC)
    N50df <- data.frame(length = n50length,
                        avgQual = n50qual,
                        GC = n50GC)
    row.names(N50df) <- "N50"
    
    # get summary
    sum <- do.call(cbind, lapply(data[, 2:4], summary))
    sum <- as.data.frame(sum, stringsAsFactors = FALSE)
    sum <- rbind(sum, N50df)
    
    # return formatted summary
    sum$avgQual <- round(as.numeric(sum$avgQual), 2)
    sum$length <- floor(as.numeric(sum$length))
    sum$GC <- round(as.numeric(sum$GC), 2)
    
    return(sum)
  })
  
  # empty plots
  pt0 <- reactive({
    if (is.null(filtered.data()))
      return(NULL)
    ggplot() + theme_void()
  })
  
  # metadata tables
  meta <- reactive({
    if (is.null(filtered.data()))
      return(NULL)
    info <- data.frame(x = c(FASTQFile()$fname,
                             nrow(filtered.data()),
                             date()))
    row.names(info) <- c("file", "seqs", "date")
    tt <- ttheme_minimal(
      base_size = 14,
      padding = unit(c(20, 6), "mm"),
      core = list(fg_params = list(hjust = 1, x = 0.9)),
      rowhead = list(fg_params = list(hjust = 1, x =
                                        0.9))
    )
    tableGrob(info, theme = tt, cols = NULL)
  })
  
  # summary table
  tbl <- reactive({
    if (is.null(filtered.summary()))
      return(NULL)
    tt <- ttheme_default(colhead = list(fg_params = list(parse = TRUE)))
    tableGrob(filtered.summary(), theme = tt)
  })
  
  # lengths
  pt1 <- reactive({
    if (is.null(filtered.data()))
      return(NULL)
    pdata <- filtered.data()
    
    # limits
    lenlims <- c(input$minlen, input$maxlen) #=c(0,1000000)
    n50length <- Nvalue(50, pdata$length)
    
    ggplot(pdata, aes(x = length)) +
      geom_density(alpha = 0.2,
                   fill = "#FF6666",
                   size = 0.2) +
      geom_vline(
        aes(xintercept = n50length),
        color = "gold3",
        linetype = "twodash",
        size = 1
      ) +
      xlim(lenlims) +
      xlab(paste(
        "read length distribution (cut at ",
        lenlims[2] / 1000,
        "kb)",
        sep = ""
      )) +
      theme_bw() +
      theme(plot.title = element_text(margin = margin(b = 0), size = 14)) +
      ggtitle(
        paste0(
          "Read length ([",
          input$minlen,
          "..",
          input$maxlen,
          "], N50=" ,
          n50length,
          ")"
        )
      )
  })
  
  # GC%
  pt2 <- reactive({
    if (is.null(filtered.data()))
      return(NULL)
    pdata <- filtered.data()
    
    # limits
    gclims <- c(input$minGC, input$maxGC)
    n50GC <- Nvalue(50, pdata$GC)
    
    ggplot(pdata, aes(x = GC)) +
      geom_density(alpha = 0.2,
                   fill = "#FF6666",
                   size = 0.2) +
      geom_vline(
        aes(xintercept = n50GC),
        color = "gold3",
        linetype = "twodash",
        size = 1
      ) +
      xlim(gclims) +
      xlab("read GC%") +
      theme_bw() +
      theme(plot.title = element_text(margin = margin(b = 0), size = 14)) +
      ggtitle(paste0("Read GC% (N50=" , n50GC, ")"))
  })
  
  # avgQual
  pt3 <- reactive({
    if (is.null(filtered.data()))
      return(NULL)
    pdata <- filtered.data()
    
    # limits
    qualims <- c(input$minqual, input$maxqual) #=c(0,40)
    n50qual <- Nvalue(50, pdata$avgQual)
    
    ggplot(pdata, aes(x = avgQual)) +
      geom_density(alpha = 0.2,
                   fill = "#FF6666",
                   size = 0.2) +
      geom_vline(
        aes(xintercept = n50qual),
        color = "gold3",
        linetype = "twodash",
        size = 1
      ) +
      xlim(qualims) +
      xlab(paste("read mean basecall qualities (cut at ", qualims[2], ")", sep =
                   "")) +
      theme_bw() +
      theme(plot.title = element_text(margin = margin(b = 0), size = 14)) +
      ggtitle(
        paste0(
          "Read BaseCall quality ([",
          input$minqual,
          "..",
          input$maxqual,
          "], N50=" ,
          n50qual,
          ")"
        )
      )
  })
  
  # bi-plot len GC%
  pt4 <- reactive({
    if (is.null(filtered.data()))
      return(NULL)
    pdata <- filtered.data()
    
    # limits
    lenlims <- c(input$minlen, input$maxlen) #=c(0,1000000)
    gclims <- c(input$minGC, input$maxGC)
    n50length <- Nvalue(50, pdata$length)
    n50GC <- Nvalue(50, pdata$GC)
    
    ggplot(pdata, aes(x = length, y = GC)) +
      geom_hex(bins = 100) +
      scale_fill_continuous(type = "viridis") +
      geom_hline(
        aes(yintercept = n50GC),
        color = "gold3",
        linetype = "twodash",
        size = 0.5
      ) +
      geom_vline(
        aes(xintercept = n50length),
        color = "gold3",
        linetype = "twodash",
        size = 0.5
      ) +
      xlim(lenlims) +
      ylim(gclims) +
      theme_bw() +
      theme(plot.title = element_text(margin = margin(b = 0), size = 14),
            legend.position = "none") +
      ggtitle("GC vs length")
  })
  
  # bi-plot len avgQual
  pt5 <- reactive({
    if (is.null(filtered.data()))
      return(NULL)
    pdata <- filtered.data()
    
    # limits
    lenlims <- c(input$minlen, input$maxlen) #=c(0,1000000)
    qualims <- c(input$minqual, input$maxqual) #=c(0,40)
    n50length <- Nvalue(50, pdata$length)
    n50qual <- Nvalue(50, pdata$avgQual)
    
    ggplot(pdata, aes(x = length, y = avgQual)) +
      geom_hex(bins = 100) +
      scale_fill_continuous(type = "viridis") +
      geom_hline(
        aes(yintercept = n50qual),
        color = "gold3",
        linetype = "twodash",
        size = 0.5
      ) +
      geom_vline(
        aes(xintercept = n50length),
        color = "gold3",
        linetype = "twodash",
        size = 0.5
      ) +
      xlim(lenlims) +
      ylim(qualims) +
      theme_bw() +
      theme(plot.title = element_text(margin = margin(b = 0), size = 14),
            legend.position = "none") +
      ggtitle("avgQual vs length")
  })
  
  # bi-plot avgQual GC%
  pt6 <- reactive({
    if (is.null(filtered.data()))
      return(NULL)
    pdata <- filtered.data()
    
    # limits
    qualims <- c(input$minqual, input$maxqual) #=c(0,40)
    gclims <- c(input$minGC, input$maxGC)
    n50qual <- Nvalue(50, pdata$avgQual)
    n50GC <- Nvalue(50, pdata$GC)
    
    ggplot(pdata, aes(x = avgQual, y = GC)) +
      geom_hex(bins = 100) +
      scale_fill_continuous(type = "viridis") +
      geom_hline(
        aes(yintercept = n50GC),
        color = "gold3",
        linetype = "twodash",
        size = 0.5
      ) +
      geom_vline(
        aes(xintercept = n50qual),
        color = "gold3",
        linetype = "twodash",
        size = 0.5
      ) +
      xlim(qualims) +
      ylim(gclims) +
      theme_bw() +
      theme(plot.title = element_text(margin = margin(b = 0), size = 14),
            legend.position = "none") +
      ggtitle("GC vs avgQual")
  })
  
  # multiplot
  plotData <- reactive({
    ptlist <-
      list(pt1(), pt2(), pt3(), tbl(), pt4(), pt5(), pt6(), meta())
    grid.arrange(grobs = ptlist,
                 nrow = 4,
                 ncol = 2)
  })
  
  output$plots <- renderPlot({
    print(plotData())
  },
  height = 1600,
  width = 800)
  
  observeEvent(input$filterFASTQ, {
    ## filter reads based on user mimits
    fun <- function(fq) {
      ## filter reads
      fq <- fq[width(fq) >= input$minlen]
      fq <- fq[width(fq) <= input$maxlen]
      fq <- fq[(alphabetScore(fq) / width(fq)) >= input$minqual]
      fq <- fq[(alphabetScore(fq) / width(fq)) <= input$maxqual]
      return(fq)
    }
    
    infq <- file.path("Uploads", "data.fq.gz")
    outfq <- file.path("Uploads", "filtered_data.fq.gz")
    fres <- filterFastq(infq, outfq, filter = fun, compress = TRUE)
    
    fastq_filter_results <- list()
    
    # filtering failed
    if (fres != outfq) {
      fastq_filter_results <-
        c(fastq_filter_results, "FAILED! filtering the FastQ data")
    }
    # all went well
    if (length(fastq_filter_results) == 0) {
      fastq_filter_results <-
        c(fastq_filter_results,
          "The filtered FastQ file was produced")
    }
    
    # report to user
    showModal(
      modalDialog(
        title = "FastQ filtering status",
        paste(unlist(fastq_filter_results), collapse = ", "),
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  
  output$downloadFASTQ <- downloadHandler(
    filename = function() {
      gsub(".fq.gz", "_filtered.fq.gz", FASTQFile()$fname)
    },
    content = function(file) {
      serverfile <- file.path("Uploads", "filtered_data.fq.gz")
      file.copy(serverfile, file)
    }
  )
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      gsub("fq.gz", input$format, FASTQFile()$fname)
    },
    content = function(file) {
      ggsave(
        file,
        plotData(),
        device = input$format,
        width = 20,
        height = 28,
        units = "cm"
      )
    }
  )
  
}

# Run the application
shinyApp(ui = ui, server = server)

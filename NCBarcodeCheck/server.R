# check two barcode files for collision
# Stephane Plaisance - VIB-NC Sep-03-2019 v0.0.1

library(shiny)
library(DNABarcodes)
library(DNABarcodeCompatibility)
library(stringr)

# load formatting functions
source("html_output.R")

###############
## functions ##
###############

define_distance_metric <- function(metric) {
  if (metric == "none") {
    return(NULL)
  }
  else{
    return(metric)
  }
}

update_valid_data <- function(df) {
  valid = c()
  if (!is.null(df)) {
    valid <- c(valid, "TRUE")
    df <- cbind(df, valid)
    return(df)
  }
  else{
    return(NULL)
  }
}

######################
## shared variables ##
######################

# distancesVec <-
#   c(
#     "Hamming" = "hamming",
#     "Seqlev" = "seqlev",
#     "Phaseshift" = "phaseshift",
#     "Levenshtein" = "levenshtein",
#     "No distance" = "none"
#   )

####################################################
## Define server logic                            ##
####################################################

shinyServer(function(input, output, session) {
  ####################################################
  ### single_search                                ###
  ####################################################
  
  ## reactive objects
  single_distance_metric <-
    reactive(define_distance_metric(input$single_distance_metric))
  single_min_distance <- reactive(input$single_min_dist)
  single_sub_cost <- reactive(input$single_sub_cost)
  single_indel_cost <- reactive(input$single_indel_cost)
  data0 <-
    reactive(update_valid_data(
      DNABarcodeCompatibility:::file_loading_and_checking(input$single_file$datapath)
    ))
  
  ### Handle when wrong input for min_distance parameter ###
  observe({
    single_min_distance <- as.numeric(single_min_distance())
    if (single_min_distance < 0 || is.na(single_min_distance)) {
      output$single_min_distance_error_message <-
        renderText("Please enter a correct positive value")
    }
    else{
      output$single_min_distance_error_message <- renderText("")
    }
  })
  
  ### Handle when wrong input for sub_cost parameter ###
  observe({
    single_sub_cost <- as.numeric(single_sub_cost())
    if (single_sub_cost < 0 || is.na(single_sub_cost)) {
      output$single_sub_cost_error_message <-
        renderText("Please enter a correct positive value")
    }
    else{
      output$single_sub_cost_error_message <- renderText("")
    }
  })
  
  ### Handle when wrong input for indel_cost parameter ###
  observe({
    single_indel_cost <- as.numeric(single_indel_cost())
    if (single_indel_cost < 0 || is.na(single_indel_cost)) {
      output$single_indel_cost_error_message <-
        renderText("Please enter a correct positive value")
    }
    else{
      output$single_min_distance_error_message <- renderText("")
    }
  })
  
  ### Disable distance numeric input when "No distance" is chosen from distance metric ###
  observeEvent(input$single_distance_metric, {
    if (is.null(single_distance_metric())) {
      disable("single_min_distance")
      updateNumericInput(session, "single_min_distance", value = 0)
    }
    else{
      enable("single_min_distance")
      updateNumericInput(session, "single_min_distance", value = 2)
    }
  })
  
  ### Load input data when single file is uploaded ###
  observeEvent(input$single_file, {
    if (is.null(data0())) {
      output$single_input_error_message <-
        renderText(paste("Error loading file : ", error_message))
      output$single_loaded_data <- renderDataTable(NULL)
    }
    else{
      output$single_input_error_message <-
        renderText(input$single_file$name)
      output$single_loaded_data <-
        renderDataTable(datatable(
          data0(),
          options = list(pageLength = 10),
          selection = "none"
        ))
      if (ncol(data0())>1) shinyjs::show("single_search")
    }
  })
  
  ####################################################
  ### When clicking on single_search button ###
  ####################################################
  
  observeEvent(input$single_search, {
    # require data in the table to avoid failed upload 
    req(data0())
    
    #Reset outputs
    output$single_table_result <- renderDataTable({
      NULL
    })
    hide("single_download_results")
    
    # set data from loaded and params
    df <- data0()

        my.dist <- function(x, y) {
      suppressWarnings(
        DNABarcodes::distance(
          x,
          y,
          metric = single_distance_metric(),
          cost_sub = single_sub_cost(),
          cost_indel = single_indel_cost()
        )
      )
    }
    #
    my.names <- function(x, y) {
      paste(x, ' <=> ', y, sep = "")
    }
    #
    
    ###############################
    # compute distances within setA
    ###############################
    lablmat <- outer(df[, 1], df[, 1], FUN = 'my.names')
    distmat <- outer(df[, 2], df[, 2], FUN = Vectorize(my.dist))
    
    # find hits
    hits <- which(distmat < single_min_distance(), arr.ind = T)
    hits <- hits[hits[, 1] != hits[, 2] , ]
    
    # find and remove mirror rows (self-pairwise comparison)
    hits <-
      hits[!duplicated(t(apply(as.data.frame(hits), 1, sort))), ]
    
    # create result table
    results <<- data.frame(
      pairs = lablmat[hits],
      distance = distmat[hits],
      setA = df[hits[, 1], 2],
      setB = df[hits[, 2], 2]
    )
    dist.label <- paste(ifelse(is.null(single_distance_metric()),"no",single_distance_metric()), "-dist", sep = "")
    colnames(results) <-
      c(
        "similar_pairs",
        dist.label,
        "barcode-A",
        "barcode-B"
      )
    
    if (!is.null(results)) {
      shinyjs::show("download_single_results")

      output$single_table_result <- renderDataTable({
        results
      })
      
      # activate  tab
      updateTabsetPanel(session, "single_tabset", selected = "single_table_result")
      }
  })
  
  ### When clicking on single_download_results button ###
  output$download_single_results <- downloadHandler(
    filename = function() {
      paste("BarcodeCollisionCheck_single_results", ".csv", sep = "")
    },
    content = function(file) {
      cat(
        paste(
          "# Pairs with small distance in the barcode set:",
          input$single_file$name,
          "\n",
          sep = " "
        ),
        file = file
      )
      write.table(
        results,
        file = file,
        row.names = FALSE,
        sep = ",",
        append = TRUE
      )
    }
  )
  
  ####################################################
  ### dual_search                                  ###
  ####################################################
  
  ## reactive objects
  distance_metric <-
    reactive(define_distance_metric(input$distance_metric))
  min_distance <- reactive(input$min_dist)
  sub_cost <- reactive(input$sub_cost)
  indel_cost <- reactive(input$indel_cost)
  data1 <-
    reactive(update_valid_data(
      DNABarcodeCompatibility:::file_loading_and_checking(input$index_file1$datapath)
    ))
  data2 <-
    reactive(update_valid_data(
      DNABarcodeCompatibility:::file_loading_and_checking(input$index_file2$datapath)
    ))
  
  ### Handle when wrong input for min_distance parameter ###
  observe({
    min_distance <- as.numeric(min_distance())
    if (min_distance < 0 || is.na(min_distance)) {
      output$min_distance_error_message <-
        renderText("Please enter a correct positive value")
    }
    else{
      output$min_distance_error_message <- renderText("")
    }
  })
  
  ### Handle when wrong input for sub_cost parameter ###
  observe({
    sub_cost <- as.numeric(sub_cost())
    if (sub_cost < 0 || is.na(sub_cost)) {
      output$sub_cost_error_message <-
        renderText("Please enter a correct positive value")
    }
    else{
      output$sub_cost_error_message <- renderText("")
    }
  })
  
  ### Handle when wrong input for indel_cost parameter ###
  observe({
    indel_cost <- as.numeric(indel_cost())
    if (indel_cost < 0 || is.na(indel_cost)) {
      output$indel_cost_error_message <-
        renderText("Please enter a correct positive value")
    }
    else{
      output$min_distance_error_message <- renderText("")
    }
  })
  
  ### Disable distance numeric input when "No distance" is chosen from distance metric ###
  observeEvent(input$distance_metric, {
    if (is.null(distance_metric())) {
      disable("min_distance")
      updateNumericInput(session, "min_distance", value = 0)
    }
    else{
      enable("min_distance")
      updateNumericInput(session, "min_distance", value = 2)
    }
  })
  
  
  ### Load input data when dual file 1 is uploaded ###
  observeEvent(input$index_file1, {
    if (is.null(data1())) {
      output$dual_input_error_message1 <-
        renderText(paste("Error loading file : ", error_message))
      output$dual_loaded_data1 <- renderDataTable(NULL)
    }
    else{
      output$dual_input_error_message1 <-
        renderText(input$index_file1$name)
      output$dual_loaded_data1 <-
        renderDataTable(datatable(
          data1(),
          options = list(pageLength = 10),
          selection = "none"
        ))
    }
    # activate tab
    updateTabsetPanel(session, "input_tab", selected = "file1")
  })
  
  ### Load input data when dual file 2 is uploaded ###
  observeEvent(input$index_file2, {
    if (is.null(data2())) {
      output$dual_input_error_message2 <-
        renderText(paste("Error loading file : ", error_message))
      output$dual_loaded_data2 <- renderDataTable(NULL)
    }
    else{
      output$dual_input_error_message2 <-
        renderText(input$index_file2$name)
      output$dual_loaded_data2 <-
        renderDataTable(datatable(
          data2(),
          options = list(pageLength = 10),
          selection = "none"
        ))
    }
    # activate tab
    updateTabsetPanel(session, "input_tab", selected = "file2")
  })
  
  ####################################################
  ### When clicking on compatibility_search button ###
  ####################################################
  observeEvent(input$compatibility_search, {
    # require data in the table to avoid failed uploads 
    req(data1(), data2())
    
    #Reset outputs
    output$dual_table_result <- renderDataTable({
      NULL
    })
    hide("dual_download_results")
    
    # set data from loaded and params
    df1 <- data1()
    df2 <- data2()
    
    my.dist <- function(x, y) {
      suppressWarnings(
        DNABarcodes::distance(
          x,
          y,
          metric = distance_metric(),
          cost_sub = sub_cost(),
          cost_indel = indel_cost()
        )
      )
    }
    #
    my.names <- function(x, y) {
      paste(x, ' <=> ', y, sep = "")
    }
    #
    
    ###############################
    # compute distances within setA
    ###############################
    lablmatA <- outer(df1[, 1], df1[, 1], FUN = 'my.names')
    distmatA <-
      outer(df1[, 2], df1[, 2], FUN = Vectorize(my.dist))
    
    # find hits
    hitsA <- which(distmatA < min_distance(), arr.ind = T)
    hitsA <- hitsA[hitsA[, 1] != hitsA[, 2] , ]
    
    # find and remove mirror rows (self-pairwise comparison)
    hitsA <-
      hitsA[!duplicated(t(apply(as.data.frame(hitsA), 1, sort))), ]
    
    # create result table (globally stored)
    resultsA <<- data.frame(
      pairs = lablmatA[hitsA],
      distance = distmatA[hitsA],
      setA = df1[hitsA[, 1], 2],
      setA = df1[hitsA[, 2], 2]
    )
    
    dist.label <- paste(ifelse(is.null(distance_metric()),"no",distance_metric()), "-dist", sep = "")
    colnames(resultsA) <-
      c(
        "similar_pairs",
        dist.label,
        "barcode-A",
        "barcode-B"
      )

    if (!is.null(resultsA)) {
      output$set1_table_result <- renderDataTable({
        resultsA
      })
    }
    
    ###############################
    # compute distances within setB
    ###############################
    lablmatB <- outer(df2[, 1], df2[, 1], FUN = 'my.names')
    distmatB <-
      outer(df2[, 2], df2[, 2], FUN = Vectorize(my.dist))
    
    # find hits
    hitsB <- which(distmatB < min_distance(), arr.ind = T)
    hitsB <- hitsB[hitsB[, 1] != hitsB[, 2] , ]
    
    # find and remove mirror rows (self-pairwise comparison)
    hitsB <-
      hitsB[!duplicated(t(apply(as.data.frame(hitsB), 1, sort))), ]
    
    # create result table (globally stored)
    resultsB <<- data.frame(
      pairs = lablmatB[hitsB],
      distance = distmatB[hitsB],
      setB = df2[hitsB[, 1], 2],
      setB = df2[hitsB[, 2], 2]
    )
    dist.label <- paste(ifelse(is.null(distance_metric()),"no",distance_metric()), "-dist", sep = "")
    colnames(resultsB) <-
      c(
        "similar_pairs",
        dist.label,
        "barcode-A",
        "barcode-B"
      )
    
    if (!is.null(resultsB)) {
      output$set2_table_result <- renderDataTable({
        resultsB
      })
    }
    ###############################
    # compute distance between sets
    ###############################
    lablmatAB <- outer(df1[, 1], df2[, 1], FUN = 'my.names')
    distmatAB <-
      outer(df1[, 2], df2[, 2], FUN = Vectorize(my.dist))
    #
    # # find hits
    hitsAB <- which(distmatAB < min_distance(), arr.ind = T)
    #
    # # create result table (globally stored)
    resultsAB <<- data.frame(
      pairs = lablmatAB[hitsAB],
      distance = distmatAB[hitsAB],
      setA = df1[hitsAB[, 1], 2],
      setB = df2[hitsAB[, 2], 2]
    )
    dist.label <- paste(ifelse(is.null(distance_metric()),"no",distance_metric()), "-dist", sep = "")
    colnames(resultsAB) <-
      c(
        "similar_pairs",
        dist.label,
        "barcode-A",
        "barcode-B"
      )
    
    if (!is.null(resultsAB)) {
      output$dual_table_result <- renderDataTable({
        resultsAB
      })
      
      shinyjs::show("download_dual_results")
      
      # activate third tab
      updateTabsetPanel(session, "dual_tabset", selected = "dual_table_result")
    }
    
  })
  
  ### When clicking on download results button ###
  output$download_dual_results <- downloadHandler(
    filename = function() {
      paste("BarcodeCollisionCheck_result", ".csv", sep = "")
    },
    content = function(file) {
      cat(
        paste(
          "# Pairs with small distance within set A:",
          input$index_file1$name,
          "\n",
          sep = " "
        ),
        file = file
      )
      write.table(
        resultsA,
        file = file,
        row.names = FALSE,
        sep = ",",
        append = TRUE
      )
      cat(
        paste(
          "# Pairs with small distance within set B:",
          input$index_file2$name,
          "\n",
          sep = " "
        ),
        file = file,
        append = TRUE
      )
      write.table(
        resultsB,
        file = file,
        row.names = FALSE,
        sep = ",",
        append = TRUE
      )
      cat(
        "# Pairs with small distance between set:",
        input$index_file1$name,
        "and set:",
        input$index_file2$name,
        "\n",
        file = file,
        append = TRUE
      )
      write.table(
        resultsAB,
        file = file,
        row.names = FALSE,
        sep = ",",
        append = TRUE
      )
    }
  )
  
  
  # the end
})

# last edit 2019-09-03 SP@NC

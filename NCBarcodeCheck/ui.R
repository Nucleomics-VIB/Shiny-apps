library(shiny)
library(shinydashboard)
library(shinyBS)
library(DT)
library(shinyjs)

sidemenuWidth <- 300

sidebar <- dashboardSidebar(width = sidemenuWidth,
                            sidebarMenu(
                                menuItem(
                                    text = "Home",
                                    tabName = "home",
                                    icon = icon("home")
                                ),
                                menuItem(
                                    text = "Tutorial",
                                    tabName = "tuto",
                                    icon = icon("graduation-cap")
                                ),
                                menuItem(
                                    text = "Single Analyse", 
                                    tabName = "single",
                                    icon=icon("wrench"),
                                    selected = TRUE
                                ),
                                menuItem(
                                    text = "Dual Analyse", 
                                    tabName = "dual",
                                    icon=icon("wrench")
                                )
                            )
                        )

body <- dashboardBody(tabItems(
    ##################
    #### TAB HOME ####
    ##################
    
    tabItem(
        tabName = "home",
        h2(
            "Barcode Collision Check"
        ),
        h3("Description"),
        p("The current Shiny-web is testing barcodes used in one or two libraries and identifies pairs that are identical or very close within and between set(s)."),
        p("Such conflicting barcodes may lead to confusion between samples during demultiplexing and require manual re-analysis of the data and read loss."),
        br(),
        h4("Status"),
        p("VIB Nucleomics Core; (Version v0.0.1, 2019-09-04)."),
        p("contributors: Stephane Plaisance,")
    ),
    
    ######################
    #### TAB TUTORIAL ####
    ######################
    
    tabItem(
        tabName = "tuto",
        h2("Tutorial"),
        tags$img(
            src = "pictureShiny.png",
            height = "70%",
            width = "70%"
        ),
        h4("1)	Upload one or two barcode tables (depending on teh chosen tool)"),
        h4("->	Each having a name column followed by a sequence column (i7)"),
        h4("2)	Set the metric and costs"),
        h4("3)	Set the minimum required distance between any barcode pair"),
        h4("PRESS Search"),
        h4("4)	Faulty barcodes in both files are returned with their distance"),
        h4("->	Columns can be sorted by clicking on their header"),
        br(),
        h4("example data 1"),
        h4("A01	CGGTTCAA"),
        h4("A02	GCTGGATT"),
        h4("A03	TAACTCGG"),
        h4("A04	TAACAGTT"),
        h4("A05	ATACTCAA"),
        h4("A06	GCTGAGAA"),
        h4("A07	ATTGGAGG"),
        h4("A08	TAGTCTAA"),
        br(),
        h4("Note: dual analysis has less distance metrics because some of them require equal number of rows and columns")
    ),
    
    #############################
    ####      PARAMETERS     ####
    #############################
    
    ####################################################
    ## single   
    ####################################################
    
    tabItem(
        tabName = "single",
        h2("Barcode Collision Check within a single pool"),
        fluidRow(
            box(
                title = "Settings",
                status = "primary",
                solidHeader = TRUE,
                width = 3,
                height = "550px",
                tipify(fileInput(
                    inputId = "single_file",
                    label = "Input file",
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"), 
                    "A file with barcode_ID and barcode_Sequence columns",
                    placement = "top"
                ),
                tipify(selectInput(
                    inputId = "single_distance_metric",
                    label = "distance metric",
                    choices = c(
                        "No distance" = "none",
                        "Hamming" = "hamming",
                        "Seqlev" = "seqlev",
                        "Phaseshift" = "phaseshift",
                        "Levenshtein" = "levenshtein"
                    ),
                    selected = "none"),
                    "method for computing the barcode distance",
                    placement = "top"
                ),
                tipify(numericInput(
                    inputId = "single_min_dist",
                    label = "min distance",
                    value = 2,
                    min = 0),
                    "pairs with distance less than this value will be shown",
                    placement = "top"
                ),
                conditionalPanel(
                    "!is.numeric(input.single_min_dist) || input.single_min_dist<0",
                    span(
                        textOutput(outputId = "single_min_distance_error_message"),
                        style = "color:red"
                    )
                ),
                numericInput(
                    inputId = "single_sub_cost",
                    label = "substitution cost",
                    value = 1,
                    min = 0
                ),
                conditionalPanel(
                    "!is.numeric(input.single_sub_cost) || input.single_sub_cost<0",
                    span(
                        textOutput(outputId = "single_sub_cost_error_message"),
                        style = "color:red"
                    )
                ),
                tipify(numericInput(
                    inputId = "single_indel_cost",
                    label = "indel cost",
                    value = 4,
                    min = 0),
                    "Set to 4 to prevent them in standard search, modify when needed",
                    placement = "top"
                ),
                conditionalPanel(
                    "!is.numeric(input.single_indel_cost) || input.single_indel_cost<0",
                    span(
                        textOutput(outputId = "single_indel_cost_error_message"),
                        style = "color:red"
                    )
                ),
                tipify(actionButton("single_search", icon = icon("search"), "Search"),
                       "click to Analyze or after changing a parameter"
                )
            ),
            tabBox(
                title = "Results: barcodes pairs with distance below threshold",
                width = 6,
                # The id lets us use input$tabset1 on the server to find the current tab
                id = "single_table_result",
                height = "350px",
                tabPanel(
                    "Within barcode Set",
                    value = "single_table_result",
                    span(
                        textOutput(outputId = "single_table_result_error_message"),
                        style = "color:red"
                    ),
                    dataTableOutput(outputId = "single_table_result"),
                    hidden(downloadButton("download_single_results", "Download results"))
                )
            )
        ),
        fluidRow(
            tabBox(
                title = "Input data",
                width = 6,
                tabPanel(
                    title = "File",
                    value = "file",
                    span(textOutput(outputId = "single_input_error_message"), style = "color:red"),
                    dataTableOutput(outputId = "single_loaded_data")
                )
            )
        )
    ),
    
    ####################################################
    ## dual
    ####################################################

    tabItem(
        tabName = "dual",
        h2("Barcode Collision Check for two pools"),
        fluidRow(
            box(
                title = "Settings",
                status = "primary",
                solidHeader = TRUE,
                width = 3,
                tipify(fileInput(
                    inputId = "index_file1",
                    label = "Input file 1",
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"), 
                    "A file with barcode_ID and barcode_Sequence columns",
                    placement = "top"
                    ),
                tipify(fileInput(
                    inputId = "index_file2",
                    label = "Input file 2",
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"), 
                    "A file with barcode_ID and barcode_Sequence columns",
                    placement = "top"
                ),
                tipify(selectInput(
                    inputId = "distance_metric",
                    label = "distance metric",
                    choices = c(
                        "No distance" = "none",
                        "Seqlev" = "seqlev"
                    ),
                    selected = "none"),
                   "method for computing the barcode distance",
                   placement = "top"
                ),
                tipify(numericInput(
                    inputId = "min_dist",
                    label = "min distance",
                    value = 2,
                    min = 0),
                    "pairs with distance less than this value will be shown",
                    placement = "top"
                ),
                conditionalPanel(
                    "!is.numeric(input.min_dist) || input.min_dist<0",
                    span(
                        textOutput(outputId = "min_distance_error_message"),
                        style = "color:red"
                    )
                ),
                numericInput(
                    inputId = "sub_cost",
                    label = "substitution cost",
                    value = 1,
                    min = 0
                ),
                conditionalPanel(
                    "!is.numeric(input.sub_cost) || input.sub_cost<0",
                    span(
                        textOutput(outputId = "sub_cost_error_message"),
                        style = "color:red"
                    )
                ),
                tipify(numericInput(
                    inputId = "indel_cost",
                    label = "indel cost",
                    value = 4,
                    min = 0),
                    "Set to 4 to prevent them in standard search, modify when needed",
                    placement = "top"
                ),
                conditionalPanel(
                    "!is.numeric(input.indel_cost) || input.indel_cost<0",
                    span(
                        textOutput(outputId = "indel_cost_error_message"),
                        style = "color:red"
                    )
                ),
                tipify(actionButton("compatibility_search", icon = icon("search"), "Search"),
                "click to Analyze or after changing a parameter"
                )
            ),
            tabBox(
                title = "Results: barcodes pairs with distance below threshold",
                width = 6,
                # The id lets us use input$tabset1 on the server to find the current tab
                id = "dual_tabset",
                height = "250px",
                tabPanel(
                    "Within Set 1",
                    value = "set1_table_result",
                    span(
                        textOutput(outputId = "set1_table_result_error_message"),
                        style = "color:red"
                    ),
                    dataTableOutput(outputId = "set1_table_result")
                ),
                tabPanel(
                    "Within Set 2",
                    value = "set2_table_result",
                    span(
                        textOutput(outputId = "set2_table_result_error_message"),
                        style = "color:red"
                    ),
                    dataTableOutput(outputId = "set2_table_result")
                ),
                tabPanel(
                    "Between Sets",
                    value = "dual_table_result",
                    span(
                        textOutput(outputId = "dual_table_result_error_message"),
                        style = "color:red"
                    ),
                    dataTableOutput(outputId = "dual_table_result"),
                    hidden(downloadButton("download_dual_results", "Download results"))
                )
            )
        ),
        fluidRow(
            tabBox(
                id = "input_tab",
                title = "Input data",
                width = 6,
                tabPanel(
                    title = "File 1",
                    value = "file1",
                    span(textOutput(outputId = "dual_input_error_message1"), style = "color:red"),
                    dataTableOutput(outputId = "dual_loaded_data1")
                ),
                tabPanel(
                    title = "File 2",
                    value = "file2",
                    span(textOutput(outputId = "dual_input_error_message2"), style = "color:red"),
                    dataTableOutput(outputId = "dual_loaded_data2")
                )
            )
        )
    )
))
ui <- fluidPage(useShinyjs(),
                dashboardPage(
                    dashboardHeader(title = "CheckBarcodeCompatibility", titleWidth = sidemenuWidth),
                    sidebar,
                    body
                ))

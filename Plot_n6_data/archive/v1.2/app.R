# Plot_n6_data.shinyapp
# A R/shiny tool to plot from n6tech exported text files
#
# Stephane Plaisance, VIB Nucleomics Core (nucleomics@vib.be)
# visit our Git: https://github.com/Nucleomics-VIB 
# version: 2024-11-27_v1.2
# © by using this tool, you accept the licence saved under ./www/licence.pdf

# Load necessary packages
library(shiny)
library(readxl)
library(ggplot2)
library(reshape2)
library(DT)
library(dplyr)
library(plotly)

# Define UI for the app
ui <- fluidPage(
  headerPanel("n6tech Raw Data Analysis"),
  windowTitle = "n6tech Raw Data Analysis",
  tags$a(href="https://nucleomicscore.sites.vib.be/en", 
         target="_blank",
         img(src="NC_logo_full.png", 
             align = "right", 
             width="150", 
             height="58.5", 
             alt="VIB Nucleomics Core"),
         h5("nucleomics@vib.be | Version: 1.2 | ", "2024-11-27")
  ), # format(Sys.Date(), "%B %d, %Y"),
  tags$a(href="license.pdf", target="_blank", "usage licence"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload Excel File", accept = c(".xlsx")),
      textInput("threshold", "Enter Fluorescence Threshold", value = "8000"),
      radioButtons("mouseOverlay", "Mouse Overlay", choices = c("Well", "Label")),
      checkboxInput("showLegend", "Show Legend", value = FALSE), # New checkbox for showing legend
      tags$br(),
      #numericInput("summaryPlotCycle", "Cycle number for summary plot:", value = 20, min = 1, max=20, step = 1),
      numericInput("summaryPlotCycle", "Cycle number for summary plot:", 
                   value = 20, min = 1, step = 1),
      uiOutput("summaryPlotCycleTooltip"),
      sliderInput("binsize", "binsize for fluo summary plots :", min = 100, max = 1000, value = 500, step = 100),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Overlay Plot", plotlyOutput("overlayPlot", height = "800px")),
        tabPanel("Sample Table",
                 textOutput("selectedCount"), # Counter for selected rows
                 actionButton("selectAllBtn", "Select All"), # Button to select all
                 actionButton("deselectAllBtn", "Deselect All"), # Button to deselect all
                 actionButton("selectDisplayedBtn", "Select Displayed Rows"),
                 actionButton("filterPlotBtn", "Filter and Redraw Plot"), # New button for filtering
                 tags$div(
                   style = "margin-top: 10px; margin-bottom: 10px;",
                   tags$p(
                     "- Search column A (A1 to A12) by typing: ", tags$code("^[A][0-9]*$"), "-> letter A followed by 1+ numbers",
                     tags$br(),
                     "- Search row 2 by typing: ", tags$code("^[A-H]*2$"), "-> letter [A to H] followed by number 2",
                     tags$br(),
                     "- Search labels with OR: ", tags$code("Pos|Neg"), "-> Labels containing 'Pos' OR 'Neg'"
                   )
                 ),
                 DTOutput("sampleTable")),
        tabPanel("Summary Table", 
                 DTOutput("summaryTable"),
                 plotOutput("stopCycleHistogram"),
                 plotOutput("fluoC25Histogram"),
                 plotOutput("endpointfluoHistogram")),
        tabPanel("Fluorescence Table", DTOutput("fluoTable"))
      )
    )
  )
)

# Define server logic for the app
server <- function(input, output, session) {
  # Reactive expression to read and process the uploaded Excel file
  data <- reactive({
    req(input$file) # Ensure a file is uploaded
    file <- input$file$datapath # Read Excel file
    raw_data <- read_excel(file, col_names = FALSE)
    
    if (nrow(raw_data) < 5) {
      return(NULL) # Return NULL if there isn't enough data
    }
    
    # Replace absent values with NA (ensure compatibility with tibbles)
    raw_data <- raw_data %>% mutate(across(everything(), ~na_if(.x, "")))
    
    # Create sample_table from the first two rows
    sample_table <- as.data.frame(t(raw_data[2:1, -c(1,ncol(raw_data))]))
    colnames(sample_table) <- c("Well", "Label")
    
    # Create fluo_table
    fluo_table <- raw_data[3:(nrow(raw_data)-3), -c(1,ncol(raw_data))]
    colnames(fluo_table) <- as.character(raw_data[2, -c(1,ncol(raw_data))])
    
    # Convert Fluorescence values to integer before reshaping
    fluo_table[] <- lapply(fluo_table, function(x) as.integer(as.numeric(as.character(x))))
    fluo_table$Cycle <- rownames(fluo_table)
    
    # Reshape fluo_table into long format
    fluo_long <- melt(fluo_table, id.vars = "Cycle", variable.name = "Well", value.name = "Fluorescence")
    
    # Add Label to the table for later choice
    fluo_long <- merge(fluo_long, sample_table, by = "Well", all.x = TRUE)
    
    # Convert Cycle to numeric if necessary
    fluo_long$Cycle <- as.numeric(as.character(fluo_long$Cycle))
    
    # Create summary_table
    summary_table <- raw_data[(nrow(raw_data) - 2):(nrow(raw_data)), -ncol(raw_data)]
    
    # Transpose summary_table and convert it back to a data frame
    summary_table <- as.data.frame(t(summary_table))
    summary_table <- summary_table[-1,]
    colnames(summary_table) <- c("Endpoint.fluo","stop.cycle","fluo.c25")
    rownames(summary_table) <- sample_table$Well
    
    summary_table$Endpoint.fluo <- as.integer(summary_table$Endpoint.fluo)
    summary_table$stop.cycle <- as.integer(summary_table$stop.cycle)
    summary_table$fluo.c25 <- as.integer(summary_table$fluo.c25)
    
    list(sample_table = sample_table, fluo_table = fluo_table, fluo_long = fluo_long, summary_table = summary_table, raw_data = raw_data)
  })
  
  # update the tooltip
  output$summaryPlotCycleTooltip <- renderUI({
    req(data())
    max_cycle <- max(data()$fluo_long$Cycle)
    tags$span(
      style = "font-size: 0.8em; color: #666;",
      title = sprintf("Allowed range: 1 to %d", max_cycle),
      icon("info-circle"), "Hover for range info"
    )
  })
  
  # Render overlay plot with mouse-over functionality using plotly and dynamic legend visibility
  output$overlayPlot <- renderPlotly({
    req(data())
    
    fluo_threshold <- as.numeric(input$threshold)
    
    overlay_text <- input$mouseOverlay 
    if (overlay_text == "Label") {
      tooltip_text <- data()$fluo_long$Label 
      legend_labels <- data()$fluo_long$Label
    } else if (overlay_text == "Well") {
      tooltip_text <- data()$fluo_long$Well 
      legend_labels <- data()$fluo_long$Well
    } else {
      tooltip_text <- "" 
      legend_labels <- data()$fluo_long$Well
    }
    
    p <- ggplot(data()$fluo_long, aes(x = Cycle, y = Fluorescence, color = Well, text = tooltip_text)) +
      geom_line(aes(group = Well)) +
      geom_point(aes(group = Well)) +
      geom_hline(yintercept = fluo_threshold, linetype = "dashed", color = "red") +
      labs(title = "Fluorescence Values Across Cycles", x = "Cycle", y = "Fluorescence", subtitle = paste("Threshold:", fluo_threshold), color = "Legend") +
      theme_minimal()
    
    ggplotly(p, tooltip = "text") %>% layout(showlegend = input$showLegend) %>%
      style(visible = TRUE, traces = seq_along(unique(data()$fluo_long$Well)))
  })
  
  # Render sample table using DT without row names and show 25 rows
  output$sampleTable <- renderDT({
    req(data())
    
    datatable(data()$sample_table,
              rownames = FALSE,
              options = list(
                pageLength = 25,
                searchHighlight = TRUE,
                search = list(regex = TRUE, caseInsensitive = TRUE),
                columnDefs = list(list(
                  searchable = TRUE,
                  targets = 0:(ncol(data()$sample_table) - 1),
                  search = list(regex = TRUE)
                ))
              ),
              filter = list(position = "top", clear = FALSE)
    )
  })

  # adapt the max value for the summaryPlotCycle after loading the data in
  observe({
    req(data())
    max_cycle <- max(data()$fluo_long$Cycle)
    updateNumericInput(session, "summaryPlotCycle", max = max_cycle)
  })
  
  # Display count of selected rows in Sample Table
  output$selectedCount <- renderText({
    paste("Selected Rows:", length(input$sampleTable_rows_selected))
  })
  
  observeEvent(input$deselectAllBtn, {
    proxy <<- dataTableProxy('sampleTable') 
    selectRows(proxy, NULL) # Deselect all rows visually
  })
  
  observeEvent(input$selectAllBtn, {
    proxy <<- dataTableProxy('sampleTable') 
    all_rows_indices <- seq_len(nrow(data()$sample_table)) # Get indices of all rows in sample table
    selectRows(proxy, all_rows_indices) # Select all rows visually without coercion issues
  })
  
  observeEvent(input$selectDisplayedBtn, {
    proxy <- dataTableProxy('sampleTable')
    
    # Get the indices of all currently displayed rows
    displayed_rows <- input$sampleTable_rows_all
    
    # Select the displayed rows
    selectRows(proxy, displayed_rows)
  })

 observeEvent(input$filterPlotBtn, {
   selected_rows_indices <- input$sampleTable_rows_selected # Get selected row indices
   
   if (length(selected_rows_indices) > 0) {
     selected_rows <- data()$sample_table[selected_rows_indices, ] # Subset of selected rows
     
     # Filter fluo_long based on selected wells/labels only
     filtered_fluo_long <- data()$fluo_long[data()$fluo_long$Well %in% selected_rows$Well, ]
     
     if (nrow(filtered_fluo_long) > 0) {
       output$overlayPlot <<- renderPlotly({
         fluo_threshold <- as.numeric(input$threshold)
         overlay_text <- input$mouseOverlay 
         
         if (overlay_text == "Label") {
           tooltip_text <- filtered_fluo_long$Label 
         } else if (overlay_text == "Well") {
           tooltip_text <- filtered_fluo_long$Well 
         } else {
           tooltip_text <- "" 
         }
         
         p_filtered <- ggplot(filtered_fluo_long, aes(x = Cycle, y = Fluorescence, color = factor(tooltip_text), text = tooltip_text)) +
           geom_line() +
           geom_point() +
           geom_hline(yintercept = fluo_threshold, linetype = "dashed", color = "red") +
           labs(title = "Fluorescence Values Across Cycles (Filtered)", x = "Cycle", y = "Fluorescence", subtitle = paste("Threshold:", fluo_threshold)) +
           theme_minimal()
         
         ggplotly(p_filtered, tooltip = "text") %>% layout(showlegend = input$showLegend)  %>%
           layout(legend = list(title = list(text = "Legend"))) # Show or hide legend based on checkbox input
       })
       
       updateTabsetPanel(session, "tabs", selected = "Overlay Plot") # Switch to Overlay Plot tab
     }
   }
 })
  
  # Render summary table using DT and show only 10 rows
  output$summaryTable <- renderDT({
    req(data())
    
    datatable(data()$summary_table,
              options = list(pageLength = 10)) # Show only 10 rows
  })
  
  # Render fluo_table using DT and show 25 rows
  output$fluoTable <- renderDT({
    req(data())
    
    datatable(data()$fluo_table,
              options = list(pageLength = 25)) # Show 25 rows
  })
  
  # Render histogram of stop.cycle in summary table tab with specified x range and labels for stop.cycle numbers
  output$stopCycleHistogram <- renderPlot({
    req(data())
    
    ggplot(data()$summary_table, aes(x = stop.cycle)) +
      geom_histogram(binwidth=1, fill="blue", color="black") +
      scale_x_continuous(limits=c(1,35), breaks=seq(1,35, by=1)) +
      labs(title="Distribution of Stop Cycle", x="Stop Cycle", y="Frequency") +
      theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12)) 
  })
  
  # Render histogram of fluo.c25 values on the Summary Table page 
  output$fluoC25Histogram <- renderPlot({
    req(data())
    
    binsize <- input$binsize
    
    ggplot(data()$summary_table, aes(x=fluo.c25)) +
      geom_histogram(binwidth=binsize, fill="red", color="black") +
      scale_x_continuous(breaks=seq(min(data()$summary_table$fluo.c25), max(data()$summary_table$fluo.c25), by=binsize)) +
      labs(title="Distribution of Fluorescence at c=25", x="Fluorescence (c=25)", y="Frequency") +
      theme(axis.text.x=element_text(size=12, angle=45,hjust=1), axis.text.y=element_text(size=12))
  })
  
  # Render histogram of Endpoint.fluo values on the Summary Table page 
  output$endpointfluoHistogram <- renderPlot({
    req(data())
    
    binsize <- input$binsize
    
    ggplot(data()$summary_table, aes(x=Endpoint.fluo)) +
      geom_histogram(binwidth=binsize, fill="green", color="black") +
      scale_x_continuous(breaks=seq(min(data()$summary_table$Endpoint.fluo), max(data()$summary_table$Endpoint.fluo), by=binsize)) +
      labs(title="Distribution of END point Fluorescence", x="END point Fluorescence", y="Frequency") +
      theme(axis.text.x=element_text(size=12, angle=45,hjust=1), axis.text.y=element_text(size=12))
  })
}

# Run the application 
shinyApp(ui=ui, server=server)
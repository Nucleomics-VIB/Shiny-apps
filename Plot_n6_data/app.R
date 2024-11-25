# Load necessary packages
library(shiny)
library(readxl)
library(ggplot2)
library(reshape2)
library(DT)
library(dplyr)
library(plotly)  # Load plotly for interactive plots

# Define UI for the app
ui <- fluidPage(
  titlePanel("n6tech Raw Data Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload Excel File", accept = c(".xlsx")),
      textInput("threshold", "Enter Fluorescence Threshold", value = "8000"),
      radioButtons("mouseOverlay", "Mouse Overlay", choices = c("Well", "Label")),  # New radio buttons
      sliderInput("binsize", 
                  "binsize for the fluo.c25 histogram :", 
                  min = 100, 
                  max = 1000, 
                  value = 500, 
                  step = 100),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Overlay Plot", plotlyOutput("overlayPlot", height = "800px")),  # Use plotlyOutput for interactivity
        tabPanel("Sample Table", DTOutput("sampleTable")),
        tabPanel("Summary Table", 
                 DTOutput("summaryTable"),  # No options here
                 plotOutput("stopCycleHistogram"),  # Existing histogram plot
                 plotOutput("fluoC25Histogram")),  # New histogram for fluo.c25
        tabPanel("Fluorescence Table", DTOutput("fluoTable"))
        #tabPanel("Raw Data Table", DTOutput("rawDataTable"))
      )
    )
  )
)

# Define server logic for the app
server <- function(input, output) {
  
  # Reactive expression to read and process the uploaded Excel file
  data <- reactive({
    req(input$file)  # Ensure a file is uploaded
    file <- input$file$datapath
    
    # Read Excel file
    raw_data <- read_excel(file, col_names = FALSE)
    
    if (nrow(raw_data) < 5) {
      return(NULL)  # Return NULL if there isn't enough data
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
    # add Label to the table for later choice
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
    
    # Ensure stop.cycle is numeric
    summary_table$stop.cycle <- as.integer(summary_table$stop.cycle)
    
    summary_table$fluo.c25 <- as.integer(summary_table$fluo.c25)
    
    list(sample_table = sample_table, fluo_table = fluo_table, fluo_long = fluo_long, summary_table = summary_table, raw_data = raw_data)
  })
  
  # Render overlay plot with mouse-over functionality using plotly and no legend
  output$overlayPlot <- renderPlotly({
    req(data())
    
    # relevant user choice in the UI
    fluo_threshold <- as.numeric(input$threshold)
    overlay_text <- input$mouseOverlay
    
    # Create a dynamic text variable based on user selection
    if (overlay_text == "Label") {
      tooltip_text <- data()$fluo_long$Label  # Assuming 'Label' is a column in your data
    } else if (overlay_text == "Well") {
      tooltip_text <- data()$fluo_long$Well   # Assuming 'Well' is a column in your data
    } else {
      tooltip_text <- ""  # Default case if needed
    }
    
    p <- ggplot(data()$fluo_long, aes(x = Cycle, y = Fluorescence, color = factor(Well), text = tooltip_text)) +
      geom_line() +
      geom_point() +
      geom_hline(yintercept = fluo_threshold, linetype = "dashed", color = "red") +
      labs(title = "Fluorescence Values Across Cycles", 
           x = "Cycle", 
           y = "Fluorescence",
           subtitle = paste("Threshold:", fluo_threshold)) +
      theme_minimal()
    
    ggplotly(p, tooltip = "text") %>% layout(showlegend = FALSE)  # Hide the legend here
  })
  
  # Render sample table using DT without row names and show 25 rows
  output$sampleTable <- renderDT({
    req(data())
    datatable(data()$sample_table, rownames = FALSE, options = list(pageLength = 25))  # Show 25 rows
  })
  
  # Render summary table using DT and show only 10 rows
  output$summaryTable <- renderDT({
    req(data())
    datatable(data()$summary_table, options = list(pageLength = 10))  # Show only 10 rows
  })
  
  # Render fluo_table using DT and show 25 rows
  output$fluoTable <- renderDT({
    req(data())
    datatable(data()$fluo_table, options = list(pageLength = 25))  # Show 25 rows
  })
  
  # Render raw_data table using DT in the new tab
  #output$rawDataTable <- renderDT({
  #  req(data())
  #  datatable(data()$raw_data)
  #})
  
  # Render histogram of stop.cycle in summary table tab with specified x range and labels for stop.cycle numbers
  output$stopCycleHistogram <- renderPlot({
    req(data())
    
    ggplot(data()$summary_table, aes(x = stop.cycle)) +
      geom_histogram(binwidth = 1, fill="blue", color="black") +
      scale_x_continuous(limits=c(1,35), breaks=seq(1,35, by=1)) + 
      labs(title = "Distribution of Stop Cycle",
           x = "Stop Cycle",
           y = "Frequency") +
      theme(axis.text.x=element_text(size=12),   # Increase x-axis text size
            axis.text.y=element_text(size=12))   # Increase y-axis text size
  })
  
  # Render histogram of fluo.c25 values on the Summary Table page
  output$fluoC25Histogram <- renderPlot({
    req(data())
    
    binsize <- input$binsize
    
    ggplot(data()$summary_table, aes(x = fluo.c25)) +
      geom_histogram(binwidth = binsize, fill="green", color="black") + 
      scale_x_continuous(breaks=seq(min(data()$summary_table$fluo.c25), max(data()$summary_table$fluo.c25), by=binsize)) +
      labs(title = "Distribution of Fluorescence at c=25",
           x = "Fluorescence (c=25)",
           y = "Frequency") +
      theme(axis.text.x=element_text(size=12), angle=45, hjust=1) +   # Rotate x-axis text to 45 degrees and increase size
      theme(axis.text.y=element_text(size=12))   # Increase y-axis text size if needed
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
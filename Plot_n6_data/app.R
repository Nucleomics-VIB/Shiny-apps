# Plot_n6_data.shinyapp
# A R/shiny tool to plot from n6tech exported text files
#
# Stephane Plaisance, VIB Nucleomics Core (nucleomics@vib.be)
# visit our Git: https://github.com/Nucleomics-VIB 
# version: 2024-11-29_v1.3
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
      fileInput("labels", "Upload Excel Labels", accept = c(".xlsx")),
      fileInput("data", "Upload Excel Data", accept = c(".xlsx")),
      textInput("threshold", "Enter Fluorescence Threshold", value = "8000"),
      radioButtons("mean_median", "Choose Statistic for Color Scale:", 
                   choices = c("Median" = "median", "Mean" = "mean"), 
                   selected = "median"),
      numericInput("cycle_number", "Focus Cycle Number:", value = 20, min = 1, max = 35),
      radioButtons("mouseOverlay", "Mouse Overlay", choices = c("Well", "Label")),
      checkboxInput("showLegend", "Show Legend", value = FALSE),
      sliderInput("binsize", "Binsize for fluo summary plots:", min = 100, max = 1000, value = 500, step = 100),
      downloadButton("downloadData", "Download Test Data"),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Overlay Plot", plotlyOutput("overlayPlot", height = "800px")),
        tabPanel("Sample Information", 
                 textOutput("selectedCount"),
                 actionButton("selectAllBtn", "Select All"),
                 actionButton("deselectAllBtn", "Deselect All"),
                 actionButton("selectDisplayedBtn", "Select Displayed Rows"),
                 actionButton("filterPlotBtn", "Filter and Redraw Plot"),
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
                 DTOutput("sampleInfo")),
        tabPanel("Summary Table", 
                 DTOutput("summaryTable"), 
                 plotOutput("stopCycleHistogram"), 
                 plotOutput("endpointfluoHistogram"), 
                 plotOutput("userCycleHistogram")),
        tabPanel("Fluorescence Plate Layout", 
                 tags$div(
                   style = "margin-top: 10px; margin-bottom: 10px;",
                   tags$p(tags$b("Note:"),": data color is centered on the statistic_value, blue cells are less fluorescent while red celles are more.",
                   tags$br(),
                   "When the End.point was reached in a former cycle, that value is repeated here")
                   ),
                 textOutput("statisticInfo"),  
                 plotOutput("fluoPlateLayout", height = "600px")),
        tabPanel("Fluorescence Table", DTOutput("fluoTable")),
        #tabPanel("melted Fluo Table", DTOutput("fluoLong")),
        #tabPanel("n6 Data", DTOutput("n6data"))
      )
    )
  )
)

# Define server logic for the app
server <- function(input, output, session) {

  # Download handler for the zip file
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("archive", ".zip", sep = "")
    },
    content = function(file) {
      full_path <- file.path("www", "archive.zip") # Path to archive.zip
      
      if (file.exists(full_path)) {
        file.copy(full_path, file) # Copy the zip file to the temp location for download
      } else {
        stop("File not found: ", full_path)
      }
    }
  )
  
  # Reactive expression to read and process the uploaded Excel NC wells2labels file
  sample_info <- reactive({
    req(input$labels)  # Ensure a file is uploaded
    labelsfile <- input$labels$datapath  # Get file path
    sample_info <- read_excel(labelsfile, col_names = TRUE)  # Read Excel file

    if (nrow(sample_info) < 5) {
      return(NULL)  # Return NULL if there isn't enough data
    }
    
    # return sample_info()
    sample_info
  })
  
  # Reactive expression to read and process the uploaded Excel n6 data file
  n6data <- reactive({
    req(input$data)  # Ensure a data file was uploaded

    datafile <- input$data$datapath  # Get file path
    n6data <- read_excel(datafile, col_names = TRUE)  # Read Excel file
    
    if (nrow(n6data) < 5) {
      return(NULL)  # Return NULL if there isn't enough data
    }
    
    # Fill NA values with the last observation carried forward
    #n6data <- na.locf(n6data, na.rm = FALSE)
    
    # return to n6data()
    n6data
  })
  
  # dynamically limit the cycle_index value based on the uploaded data
  max_cycle <- reactive({
    req(n6data())
    nrow(n6data())-2
  })
  
  observe({
    req(max_cycle())
    updateNumericInput(session, "cycle_number", max = max_cycle())
    
    # Check if the current value exceeds the max and update if necessary
    if (input$cycle_number > max_cycle()) {
      updateNumericInput(session, "cycle_number", value = max_cycle())
    }
  })
  
  # add labels where needed when both files are uploaded
  data <- reactive({
    req(sample_info(), n6data(), input$cycle_number)
    
    # copy external objects
    n6data <- n6data()
    sample_info <- sample_info()
    
    # split n6data into fluo_data and summary_data  
    
    ###################################################
    # remove last two rows of raw data to get fluo_data
    fluo_table <- n6data[-c((nrow(n6data)-1):(nrow(n6data))),]
    
    ## Convert Fluorescence values to integer before reshaping
    fluo_table <- as.data.frame(lapply(fluo_table, function(x) as.integer(as.numeric(as.character(x)))))
    
    ## Reshape fluo_table into long format
    fluo_long <- melt(fluo_table, id.vars = "Cycle", variable.name = "Well", value.name = "Fluorescence")
    
    # add labels from sample_info
    fluo_long <- merge(fluo_long, sample_info[, c("Well", "Label")], by = "Well", all.x = TRUE)
    
    ########################################################
    # Create summary_table from the last two rows of n6data
    summary_table <- rbind(colnames(n6data), n6data[(nrow(n6data)-1):(nrow(n6data)),])
    
    # collect initial Cycle20 values
    #debugging
    #cycle_index <- 20
    cycle_index <- input$cycle_number
    row_title <- paste("Cycle", cycle_index, sep=".")
    
    # Match cycle number with fluo table and extract corresponding row values
    matched_row <- as.data.frame(fluo_table[fluo_table$Cycle == cycle_index, -1])  # Exclude 'Cycle' column
    rownames(matched_row) <- row_title
    
    # Transpose summary_table and convert it back to a data frame
    summary_table <- as.data.frame(rbind(summary_table[,-1], matched_row))
    rownames(summary_table) <- c("Well","Endpoint.fluo","stop.cycle", row_title) # Rename columns
    summary_table <- as.data.frame(t(summary_table))
    
    # convert to integers
    summary_table[, 2:4] <- lapply(summary_table[, 2:4], as.integer)
    
    # add Labels
    summary_table <- merge(summary_table, sample_info[, c("Well", "Label")], by = "Well", all.x = TRUE)
    
    # output tables readable as data()$<tablename>
    list(fluo_table = fluo_table, fluo_long = fluo_long, summary_table = summary_table)
  })
  
  # Externalized overlayPlot plotting function
  create_overlay_plot <- function(fluo_data, fluo_threshold, overlay_text, legend_group, tooltip_text, show_legend) {
    p <- ggplot(fluo_data, aes(x = Cycle, y = Fluorescence, color = factor(legend_group), text = tooltip_text)) +
      geom_line() +
      geom_point() +
      geom_hline(yintercept = fluo_threshold, linetype = "dashed", color = "red") +
      labs(title = "Fluorescence Values Across Cycles", 
           x = "Cycle", 
           y = "Fluorescence", 
           subtitle = paste("Threshold:", fluo_threshold), 
           color = overlay_text) +
      theme_minimal()
    
    ggplotly(p, tooltip = "text") %>% 
      layout(showlegend = show_legend) %>%
      layout(legend = list(title = list(text = overlay_text)))
  }
  
  # Render the original overlay plot with mouse-over functionality using plotly and dynamic legend visibility
  output$overlayPlot <- renderPlotly({
    req(data()$fluo_long)
    fluo_threshold <- as.numeric(input$threshold)
    overlay_text <- input$mouseOverlay 
    
    # Create a dynamic text variable based on user selection
    if (overlay_text == "Label") {
      tooltip_text <- data()$fluo_long$Label 
      legend_group <- data()$fluo_long$Label
    } else {
      tooltip_text <- data()$fluo_long$Well 
      legend_group <- data()$fluo_long$Well
    }
    
    create_overlay_plot(data()$fluo_long, fluo_threshold, overlay_text, legend_group, tooltip_text, input$showLegend)
  })
  
  # Render sampleInfo using DT without row names and show 25 rows
  output$sampleInfo <- renderDT({
    req(sample_info())
    
    datatable(sample_info(),
              rownames = FALSE,
              options = list(
                pageLength = 25,
                searchHighlight = TRUE,
                search = list(regex = TRUE, caseInsensitive = TRUE),
                columnDefs = list(list(
                  searchable = TRUE,
                  targets = 0:(ncol(sample_info()) - 1),
                  search = list(regex = TRUE)
                ))
              ),
              filter = list(position = "top", clear = FALSE)
    )
  })

  # Display count of selected rows in sample_info Table
  output$selectedCount <- renderText({
    paste("Selected Rows:", length(input$sampleInfo_rows_selected))
  })
  
  observeEvent(input$deselectAllBtn, {
    proxy <<- dataTableProxy('sampleInfo') 
    selectRows(proxy, NULL) # Deselect all rows visually
  })
  
  observeEvent(input$selectAllBtn, {
    proxy <<- dataTableProxy('sampleInfo') 
    all_rows_indices <- seq_len(nrow(sample_info())) # Get indices of all rows in sample table
    selectRows(proxy, all_rows_indices) # Select all rows visually without coercion issues
  })
  
  observeEvent(input$selectDisplayedBtn, {
    proxy <- dataTableProxy('sampleInfo')
    
    # Get the indices of all currently displayed rows
    displayed_rows <- input$sampleInfo_rows_all
    
    # Select the displayed rows
    selectRows(proxy, displayed_rows)
  })
  
  observeEvent(input$filterPlotBtn, {
    req(data()$fluo_long, input$sampleInfo_rows_selected)
    selected_rows_indices <- input$sampleInfo_rows_selected
    
    if (length(selected_rows_indices) > 0) {
      selected_rows <- sample_info()[selected_rows_indices, ]
      
      filtered_fluo_long <- data()$fluo_long[data()$fluo_long$Well %in% selected_rows$Well, ]
      
      if (nrow(filtered_fluo_long) > 0) {
        output$overlayPlot <- renderPlotly({
          fluo_threshold <- as.numeric(input$threshold)
          overlay_text <- input$mouseOverlay 
          
          legend_group <- if (overlay_text == "Label") filtered_fluo_long$Label else filtered_fluo_long$Well
          tooltip_text <- legend_group
          
          create_overlay_plot(filtered_fluo_long, fluo_threshold, overlay_text, legend_group, tooltip_text, input$showLegend)
        })
        
        updateTabsetPanel(session, "tabs", selected = "Overlay Plot")
      }
    }
  })
  
  # render summary_table
  output$summaryTable <- renderDT({
    req(data()$summary_table)
    datatable(data()$summary_table, rownames = FALSE)
  })
  
  # render full fluo_data
  output$fluoTable <- renderDT({
    req(data()$fluo_table)
    datatable(data()$fluo_table, rownames = FALSE, 
              options = list(pageLength = 50))
  })
  
  # render plot for stop.cycle
  output$stopCycleHistogram <- renderPlot({
    req(data()$summary_table)
    ggplot(data()$summary_table, aes(x=stop.cycle)) +
      geom_histogram(binwidth=1, fill="blue", color="black") +
      scale_x_continuous(limits=c(1,35), breaks=seq(1,35, by=1)) + 
      labs(title = "Distribution of Stop Cycle",
           x = "Stop Cycle",
           y = "Frequency") +
      theme_minimal() + 
      theme(axis.text.x=element_text(size=12),   # Increase x-axis text size
            axis.text.y=element_text(size=12))   # Increase y-axis text size
  })
  
  # render plot for Endpoint.fluo
  output$endpointfluoHistogram <- renderPlot({
    req(data()$summary_table)
    
    binsize <- input$binsize
    
    ggplot(data()$summary_table, aes(x=Endpoint.fluo)) +
      geom_histogram(binwidth=binsize, fill="green", color="black") +
      scale_x_continuous(breaks=seq(min(data()$summary_table$Endpoint.fluo), max(data()$summary_table$Endpoint.fluo), by=binsize)) +
      labs(title="Distribution of END point Fluorescence", x="END point Fluorescence", y="Frequency") +
      theme_minimal() + 
      theme(axis.text.x=element_text(size=12, angle=45,hjust=1), axis.text.y=element_text(size=12))
  })
  
  # render plot for userCycleHistogram
  output$userCycleHistogram <- renderPlot({
    req(data()$summary_table)
    cycle_index <- input$cycle_number
    column_name <- paste("Cycle", cycle_index, sep=".")
    binsize <- input$binsize
    
    # Check if the column exists in the summary_table
    if (!(column_name %in% names(data()$summary_table))) {
      return(NULL)  # Return NULL if the column doesn't exist
    }
    
    #  # Extract cycle data from the full fluo_table
    fluo_data_for_cycle <- data()$fluo_table[data()$fluo_table$Cycle == cycle_index, -1]
      
    # Replace NA values with endpoint fluorescence
    endpoint_fluo <- setNames(data()$summary_table$Endpoint.fluo, data()$summary_table$Well)
    fluo_data_for_cycle <- sapply(names(fluo_data_for_cycle), function(well) {
        value <- fluo_data_for_cycle[[well]]
        if (is.na(value)) endpoint_fluo[well] else value
      })
      
    # Create a data frame for plotting
    plot_data <- data.frame(Fluorescence = unlist(fluo_data_for_cycle))
    
    ggplot(plot_data, aes(x = Fluorescence)) +
      geom_histogram(binwidth = binsize, fill = "red", color = "black") +
      labs(title = paste0("Distribution of Fluorescence at Cycle: ", cycle_index),
           x = "Fluorescence",
           y = "Frequency") +
      theme_minimal() +
      theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 12))
  })
  
  # write the value of statistic_value in the heatmap tab
  output$statisticInfo <- renderText({
    req(data()$fluo_table)
    cycle_index <- input$cycle_number
    fluo_data_for_cycle <- data()$fluo_table[data()$fluo_table$Cycle == cycle_index, -1]
    
    statistic_type <- input$mean_median
    statistic_value <- if (statistic_type == "median") {
      median(unlist(fluo_data_for_cycle), na.rm = TRUE)
    } else {
      round(mean(unlist(fluo_data_for_cycle), na.rm = TRUE))
    }
    
    paste(toupper(statistic_type), " = ", statistic_value)
  })
  
  output$fluoPlateLayout <- renderPlot({
    req(data()$fluo_table, data()$summary_table)
    cycle_index <- input$cycle_number
    
    # Extract cycle data from the full fluo_table
    fluo_data_for_cycle <- data()$fluo_table[data()$fluo_table$Cycle == cycle_index, -1]
    
    # Create a named vector of endpoint fluorescence values
    endpoint_fluo <- setNames(data()$summary_table$Endpoint.fluo, data()$summary_table$Well)
    
    # Replace NA values with endpoint fluorescence
    fluo_data_for_cycle <- sapply(names(fluo_data_for_cycle), function(well) {
      value <- fluo_data_for_cycle[[well]]
      if (is.na(value)) endpoint_fluo[well] else value
    })
    
    heatmap_data <- data.frame(
      Row = as.factor(rep(LETTERS[1:8], each = 12)),
      Column = rep(1:12, times = 8),
      Fluorescence = unlist(fluo_data_for_cycle)
    )
    
    statistic_value <- if (input$mean_median == "median") {
      median(heatmap_data$Fluorescence, na.rm = TRUE)
    } else {
      mean(heatmap_data$Fluorescence, na.rm = TRUE)
    }
    
    min_fluo <- min(heatmap_data$Fluorescence, na.rm = TRUE)
    max_fluo <- max(heatmap_data$Fluorescence, na.rm = TRUE)
    
    p <- ggplot(heatmap_data, aes(x = Column, y = Row, fill = Fluorescence)) +
      geom_tile(color = "black") +
      scale_fill_gradientn(
        colors = c("darkblue", "white", "red"),
        values = scales::rescale(c(min_fluo, statistic_value, max_fluo)),
        limits = c(min_fluo, max_fluo),
        name = "Fluorescence"
      ) +
      labs(title = paste("Fluorescence Plate Layout at Cycle", cycle_index), y = "Rows") +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 20)
      )
    
    p + scale_x_continuous(
      breaks = 1:12,
      sec.axis = sec_axis(~., breaks = 1:12, labels = 1:12, name = NULL)
    ) +
      theme(
        axis.text.y = element_text(size = 14),
        axis.text.x.top = element_text(size = 14),
        axis.ticks.x.top = element_line(size = 0.5)
      ) +
      scale_y_discrete(limits = rev(levels(heatmap_data$Row)))
  })
  
  # render n6data (debugging)
  #output$n6data <- renderDT({
  #  req(n6data())
  #  datatable(n6data(), rownames = FALSE)
  #})

  # render data()$fluo_long (debugging)  
  #output$fluoLong <- renderDT({
  #  req(data()$fluo_long)
  #  datatable(data()$fluo_long, rownames = FALSE)
  #})
  
}

shinyApp(ui=ui, server=server)
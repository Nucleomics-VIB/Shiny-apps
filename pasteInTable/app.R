library(shiny)
library(DT)
library(dplyr)

# Replace region in a table by user-pasted new data
# SP@NC 2024-12-09

# Sample initial data.frame
initial_data <- data.frame(
  A = 1:5,
  B = 6:10,
  C = 11:15,
  D = 16:20
)

ui <- fluidPage(
  h3("Editable Table with Region Replacement (Including Empty Cells)"),
  
  # Editable table
  DTOutput("table"),
  
  # Instructions and inputs
  p("Select a cell in the table (upper-left corner of target region), paste your data below, and click 'Replace Region'."),
  p("Use '\\t' for empty cells. Empty lines will be skipped."),
  textAreaInput("pasted_data", "Paste your data here (tab-separated):", rows = 5),
  actionButton("replace_btn", "Replace Region")
)

server <- function(input, output, session) {
  
  # Reactive value to store the data.frame
  rv <- reactiveValues(data = initial_data, selected_cell = NULL)
  
  # Render editable DataTable
  output$table <- renderDT({
    datatable(rv$data, editable = TRUE, selection = list(mode = "single", target = "cell"))
  })
  
  # Observe cell selection
  observeEvent(input$table_cell_clicked, {
    info <- input$table_cell_clicked
    if (!is.null(info$row) && !is.null(info$col)) {
      rv$selected_cell <- c(row = info$row, col = info$col)
    }
  })
  
  # Replace the selected region with pasted data when button is clicked
  observeEvent(input$replace_btn, {
    req(input$pasted_data) # Ensure something is pasted
    req(rv$selected_cell) # Ensure a cell is selected
    
    # Parse pasted data (assuming tab-separated values)
    pasted_data <- strsplit(input$pasted_data, "\n")[[1]] %>%
      lapply(function(row) {
        cells <- strsplit(row, "\t")[[1]]
        cells[cells == ""] <- NA  # Convert empty strings to NA
        return(cells)
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame(stringsAsFactors = FALSE)
    
    # Remove empty rows
    pasted_data <- pasted_data[rowSums(!is.na(pasted_data)) > 0, , drop = FALSE]
    
    # Convert pasted data to numeric if possible, keeping NA values
    pasted_data[] <- lapply(pasted_data, function(x) {
      as.numeric(ifelse(is.na(x), NA, x))
    })
    
    # Get dimensions of pasted data
    nrows <- nrow(pasted_data)
    ncols <- ncol(pasted_data)
    
    # Get starting position from selected cell
    start_row <- rv$selected_cell["row"]
    start_col <- rv$selected_cell["col"]
    
    # Replace part of the existing table with pasted data
    for (i in seq_len(nrows)) {
      for (j in seq_len(ncols)) {
        if ((start_row + i - 1) <= nrow(rv$data) && (start_col + j - 1) <= ncol(rv$data)) {
          rv$data[start_row + i - 1, start_col + j - 1] <- pasted_data[i, j]
        }
      }
    }
    
    # Update the table
    rv$data <- rv$data
  })
}

shinyApp(ui, server)
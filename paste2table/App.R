library(shiny)
library(DT)
library(readr)

# Create table from user-pasted 2D data (eg. excel copy)
# SP@NC 2024-12-09

ui <- fluidPage(
  titlePanel("Excel Data Import"),
  
  sidebarLayout(
    sidebarPanel(
      textAreaInput("data_input", "Paste your Excel data here:", rows = 10),
      checkboxInput("header", "First row as header", value = TRUE),
      selectInput("sep", "Field separator:",
                  choices = c("Tab" = "\t", "Comma" = ",", "Semicolon" = ";"),
                  selected = "\t"),
      selectInput("quote", "Quote character:",
                  choices = c("None" = "", "Double quote" = '"', "Single quote" = "'"),
                  selected = '"'),
      numericInput("skip", "Skip rows:", value = 0, min = 0),
      actionButton("update", "Update Table")
    ),
    
    mainPanel(
      DTOutput("table")
    )
  )
)

server <- function(input, output, session) {
  data <- reactiveVal(NULL)
  
  observeEvent(input$update, {
    tryCatch({
      df <- read_delim(input$data_input, 
                       delim = input$sep, 
                       quote = input$quote,
                       col_names = input$header,
                       skip = input$skip)
      data(df)
    }, error = function(e) {
      showNotification("Error parsing data. Please check your input and options.", type = "error")
    })
  })
  
  output$table <- renderDT({
    req(data())
    datatable(data(), options = list(scrollX = TRUE))
  })
}

shinyApp(ui, server)
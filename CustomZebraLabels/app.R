#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # label type
    # if medium or big is selected, only row1 will be printed
    radioButtons("labeltype", label = h3("Label Type"),
                 choices = list("text" = 1, "numBarcode" = 2, "textBarcode" = 3), 
                 selected = 1),
    
    # font size
    # if medium or big is selected, only row1 will be printed
    radioButtons("fontsize", label = h3("select Font size (1 row)"),
                 choices = list("normal" = 1, "medium" = 2, "big" = 3), 
                 selected = 1),
    
    hr(),
    
    # text rows
    h3("Enter text for 1 to 5 rows"),
    textInput("text1", label="", value = "Enter text..."),
    textInput("text2", label="", value = ""),
    textInput("text3", label="", value = ""),
    textInput("text4", label="", value = ""),
    textInput("text5", label="", value = ""),

    # copies
    sliderInput("copies", label = h3("Copies"), min = 1, 
                max = 10, value = 1),
    
    # Copy the line below to make an action button
    actionButton("print", label = h2("print")),

)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    # resuls
    output$value1 <- renderPrint({ input$text1 })
    output$value2 <- renderPrint({ input$text2 })
    output$value3 <- renderPrint({ input$text3 })   
    output$value4 <- renderPrint({ input$text4 })
    output$value5 <- renderPrint({ input$text5 })
}

# Run the application 
shinyApp(ui = ui, server = server)

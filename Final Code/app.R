library(shiny)

ui <- fluidPage(
  titlePanel("Gene Expression Data Preprocessing"),
  
  sidebarLayout(
    sidebarPanel(
      # Dataset selection
      radioButtons("dataset", "Choose Dataset:",
                   choices = list("Rat Data (GDS2901)" = "rat",
                                  "Golub Data" = "golub"),
                   selected = "rat"),
      
      # File path inputs
      conditionalPanel(
        condition = "input.dataset == 'rat'",
        textInput("rat_path", "GDS2901.soft file path:", 
                  value = "data/GDS2901.soft")
      ),
      conditionalPanel(
        condition = "input.dataset == 'golub'",
        textInput("golub_path", "Golub.txt file path:", 
                  value = "data/golub.txt")
      ),
      
      # Process button
      actionButton("process", "Process Data", 
                   class = "btn-primary")
    ),
    
    mainPanel(
      # Output messages
      verbatimTextOutput("status"),
      verbatimTextOutput("log")
    )
  )
)

server <- function(input, output, session) {
  # Create reactive values to store processing status
  values <- reactiveValues(
    status = "",
    log = ""
  )
  
  # Process button handler
  observeEvent(input$process, {
    values$status <- "Processing..."
    
    # Determine which file path to use
    file_path <- if(input$dataset == "rat") {
      input$rat_path
    } else {
      input$golub_path
    }
    
    # Run preprocessing
    tryCatch({
      result <- system2("./preprocess",
                        args = c(input$dataset, file_path),
                        stdout = TRUE,
                        stderr = TRUE,
                        wait = TRUE)
      
      values$status <- "Processing complete!"
      values$log <- paste(result, collapse = "\n")
      
    }, error = function(e) {
      values$status <- "Error!"
      values$log <- paste("Error processing data:", e$message)
    })
  })
  
  # Display status and log
  output$status <- renderText({
    values$status
  })
  
  output$log <- renderText({
    values$log
  })
}

shinyApp(ui = ui, server = server)
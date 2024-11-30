library(shiny)

ui <- fluidPage(
  titlePanel("Gene Expression Data Preprocessing and Clustering"),
  
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
      
      # Process and Cluster buttons
      actionButton("process", "Process Data", 
                   class = "btn-primary"),
      br(), br(),  # Add some spacing
      actionButton("cluster", "Perform Clustering",
                   class = "btn-success")
    ),
    
    mainPanel(
      # Status and processing outputs
      h4("Processing Status:"),
      verbatimTextOutput("status"),
      verbatimTextOutput("log"),
      
      # Clustering outputs
      h4("Clustering Results:"),
      verbatimTextOutput("clustering_info"),
      
      # Plots
      h4("Clustering Visualizations:"),
      plotOutput("cluster_comparison_plot"),
      plotOutput("diffcoex_distribution_plot"),
      plotOutput("coxpress_distribution_plot"),
      plotOutput("cluster_overlap_plot")
    )
  )
)

server <- function(input, output, session) {
  # Create reactive values to store processing status
  values <- reactiveValues(
    status = "",
    log = "",
    preprocessing_complete = FALSE,
    clustering_complete = FALSE,
    clustering_plots = NULL,
    clustering_info = NULL
  )
  
  # Process button handler
  observeEvent(input$process, {
    values$status <- "Processing..."
    values$preprocessing_complete <- FALSE
    values$clustering_complete <- FALSE
    
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
      
      values$status <- "Processing complete! Ready for clustering."
      values$log <- paste(result, collapse = "\n")
      values$preprocessing_complete <- TRUE
      
    }, error = function(e) {
      values$status <- "Error!"
      values$log <- paste("Error processing data:", e$message)
      values$preprocessing_complete <- FALSE
    })
  })
  
  # Cluster button handler
  observeEvent(input$cluster, {
    values$status <- "Clustering..."
    
    tryCatch({
      # Get the appropriate file paths based on dataset
      if (input$dataset == "rat") {
        data_paths <- list(
          diffcoex = list(
            condition1 = "output/diffcoex/rat_wild_types.csv",
            condition2 = "output/diffcoex/rat_eker_mutants.csv"
          ),
          coxpress = list(
            condition1 = "output/coxpress/rat_wild_types.csv",
            condition2 = "output/coxpress/rat_eker_mutants.csv"
          )
        )
      } else {
        data_paths <- list(
          diffcoex = list(
            condition1 = "output/diffcoex/golub_ALL_samples.csv",
            condition2 = "output/diffcoex/golub_AML_samples.csv"
          ),
          coxpress = list(
            condition1 = "output/coxpress/golub_ALL_samples.csv",
            condition2 = "output/coxpress/golub_AML_samples.csv"
          )
        )
      }
      
      # Source the clustering script and perform clustering
      source("clustering.R")
      results <- performClustering(data_paths)
      
      # Store results in reactive values
      values$clustering_complete <- TRUE
      values$clustering_plots <- results$plots
      values$clustering_info <- paste(
        "DiffCoEx Summary:",
        "Total Genes:", results$diffcoex_summary$Total_Genes,
        "\nTotal Clusters:", results$diffcoex_summary$Total_Clusters,
        "\nGenes Per Cluster Summary:",
        "\n  Min:", results$diffcoex_summary$Genes_Per_Cluster[1],
        "\n  1st Qu:", results$diffcoex_summary$Genes_Per_Cluster[2],
        "\n  Median:", results$diffcoex_summary$Genes_Per_Cluster[3],
        "\n  Mean:", round(results$diffcoex_summary$Genes_Per_Cluster[4], 2),
        "\n  3rd Qu:", results$diffcoex_summary$Genes_Per_Cluster[5],
        "\n  Max:", results$diffcoex_summary$Genes_Per_Cluster[6],
        "\n\nCoXpress Summary:",
        "\nTotal Genes:", results$coxpress_summary$Total_Genes,
        "\nTotal Clusters:", results$coxpress_summary$Total_Clusters,
        "\nGenes Per Cluster Summary:",
        "\n  Min:", results$coxpress_summary$Genes_Per_Cluster[1],
        "\n  1st Qu:", results$coxpress_summary$Genes_Per_Cluster[2],
        "\n  Median:", results$coxpress_summary$Genes_Per_Cluster[3],
        "\n  Mean:", round(results$coxpress_summary$Genes_Per_Cluster[4], 2),
        "\n  3rd Qu:", results$coxpress_summary$Genes_Per_Cluster[5],
        "\n  Max:", results$coxpress_summary$Genes_Per_Cluster[6],
        sep = " "
      )
      values$status <- "Clustering complete!"
      
    }, error = function(e) {
      values$status <- "Clustering Error!"
      values$log <- paste("Error during clustering:", e$message)
      values$clustering_complete <- FALSE
    })
  })
  
  # Display status and log
  output$status <- renderText({
    values$status
  })
  
  output$log <- renderText({
    values$log
  })
  
  # Display clustering results
  output$clustering_info <- renderText({
    req(values$clustering_complete)
    values$clustering_info
  })
  
  output$cluster_comparison_plot <- renderPlot({
    req(values$clustering_complete)
    values$clustering_plots$cluster_comparison
  })
  
  output$diffcoex_distribution_plot <- renderPlot({
    req(values$clustering_complete)
    ggplot(values$clustering_plots$diffcoex_distribution, 
           aes(x = Cluster, y = Gene_Count)) +
      geom_bar(stat = "identity", fill = "blue", alpha = 0.7) +
      labs(title = "Number of Genes Per Cluster - DiffCoEx", 
           x = "Cluster", y = "Number of Genes") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            text = element_text(size = 14))
  })
  
  output$coxpress_distribution_plot <- renderPlot({
    req(values$clustering_complete)
    ggplot(values$clustering_plots$coxpress_distribution, 
           aes(x = Cluster, y = Gene_Count)) +
      geom_bar(stat = "identity", fill = "orange", alpha = 0.7) +
      labs(title = "Number of Genes Per Cluster - CoXpress", 
           x = "Cluster", y = "Number of Genes") +
      theme_minimal() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            text = element_text(size = 14))
  })
  
  output$cluster_overlap_plot <- renderPlot({
    req(values$clustering_complete)
    values$clustering_plots$cluster_overlap
  })
}

shinyApp(ui = ui, server = server)
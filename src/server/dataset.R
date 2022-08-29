load_dataset <- reactive({
  path <- file.path("..", "data", "datasets", paste(input$dataset_name, ".csv", sep=""))
  read.csv(path)
})

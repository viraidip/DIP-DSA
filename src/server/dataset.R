load_dataset <- reactive({
  path <- file.path("..", "data", paste(input$dataset_name, ".csv", sep=""))
  read.csv(path)
})

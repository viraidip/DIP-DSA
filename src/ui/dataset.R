dataset_tab <- tabItem(tabName="dataset",
  h1("Select/Load dataset"),
  fluidRow(
    box(
      title="Select dataset",
      selectInput(
        inputId="dataset_name",
        label="select an existing dataset:",
        choices=c("California2007")
      ),
    ),
    box(
      title="upload dataset",
      "Allows the user to upload a custom dataset"
    ),
    box(
      title="overview",
      width=12,
      "display a table with the data",
      tableOutput("dataset_table")
    )
  )
)

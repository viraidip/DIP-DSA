library(tools)

all_datasets <- tools::file_path_sans_ext(list.files(file.path("..", "data", "datasets")))

dataset_tab <- tabItem(tabName="dataset",
  h1("Select/Load dataset"),
  fluidRow(
    box(
      title="Select dataset",
      selectInput(
        inputId="strain",
        label="select an existing dataset:",
        choices=all_datasets
      ),
    ),
    box(
      title="upload dataset",
      "Allows the user to upload a custom dataset",
      "Not only the dataset but also the full length sequences have to be uploaded",
      fileInput(
        inputId="dataset_file",
        label="Upload custom dataset here:"
      )
    ),
    box(
      title="overview",
      width=12,
      "display a table with the data",
      dataTableOutput("dataset_table")
    )
  )
)

library(tools)

all_datasets <- tools::file_path_sans_ext(list.files(file.path("..", "data", "datasets")))

dataset_tab <- tabItem(tabName="dataset",
  h1("Select an existing dataset or load a custom one"),
  fluidRow(
    box(
      title="Select dataset",
      dropMenu(
        dropdownButton("Info", status='success', icon=icon('info')),
        h3(strong('Information')),
        br(),
        h5('This is really helpful'),
      ),
      selectInput(
        inputId="strain",
        label="select an existing dataset:",
        choices=all_datasets
      ),
    ),
    box(
      title="Upload dataset",
      "Allows the user to upload a custom dataset",
      "Not only the dataset but also the full length sequences have to be uploaded",
      fileInput(
        inputId="dataset_file",
        label="Upload custom dataset here:"
      )
    ),
    box(
      title="Overview",
      width=12,
      "display a table with the data",
      dataTableOutput("dataset_table")
    )
  )
)

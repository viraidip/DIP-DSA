library(tools)

all_datasets <- tools::file_path_sans_ext(list.files(DATASETSPATH))

dataset_tab <- tabItem(tabName="dataset",
  h1("Select an existing dataset or load a custom one"),
  fluidRow(
    box(
      title="Select dataset",
      selectInput(
        inputId="strain",
        label="select an existing dataset:",
        choices=all_datasets
      ),
      "More info about the datasets can be found in the 'About' tab."
    ),
    box(
      width=12,
      title="Upload dataset",
      textInput(
        inputId="upload_strain",
        label="Write strain name"
      ),
      fileInput(
        inputId="upload_dataset_file",
        label="Upload custom dataset here (.csv, .xlsx):",
        accept=c(".csv")
      ),
      div(style="display:inline-block",
        fileInput(
          inputId="upload_PB2_file",
          label="Upload PB2 sequence as FASTA:"
        )
      ),
      div(style="display:inline-block",
        fileInput(
          inputId="upload_PB1_file",
          label="Upload PB1 sequence as FASTA:"
        )
      ),
      div(style="display:inline-block",
        fileInput(
          inputId="upload_PA_file",
          label="Upload PA sequence as FASTA:"
        )
      ),
      div(style="display:inline-block",
        fileInput(
          inputId="upload_HA_file",
          label="Upload HA sequence as FASTA:"
        )
      ),
      div(style="display:inline-block",
        fileInput(
          inputId="upload_NP_file",
          label="Upload NP sequence as FASTA:"
        )
      ),
      div(style="display:inline-block",
        fileInput(
          inputId="upload_NA_file",
          label="Upload NA sequence as FASTA:"
        )
      ),
      div(style="display:inline-block",
        fileInput(
          inputId="upload_M_file",
          label="Upload M sequence as FASTA:"
        )
      ),
      div(style="display:inline-block",
        fileInput(
          inputId="upload_NS_file",
          label="Upload NS sequence as FASTA:"
        )
      ),
      br(),
      actionButton(
        inputId="dataset_submit",
        label="Upload data"
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

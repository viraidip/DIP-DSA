new_dataset_tab <- tabItem(tabName="new_dataset",
  h1("Enter new data"),
  fluidRow(
    box(
      width=12,
      title="Upload new dataset",
      selectizeInput(
        inputId="upload_strain",
        label="Select an existing strain:",
        choices=gsub(
          "_",
          "/",
          list.dirs(DATASETSPATH, full.names=FALSE, recursive=FALSE)
        )
      ),
      textInput(
        inputId="upload_dataset",
        label="Give the dataset a unique name:"
      ),
      fileInput(
        inputId="upload_dataset_file",
        label="Upload a custom dataset (*.csv):",
        accept=c(".csv")
      ),
      br(),
      actionButton(
        inputId="dataset_submit",
        label="Upload data"
      )
    ),
    box(
      width=12,
      title="Create new strain",
      textInput(
        inputId="new_strain",
        label="Write the name of the new strain:",
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
        inputId="strain_submit",
        label="Create new strain"
      )
    )
  )
)

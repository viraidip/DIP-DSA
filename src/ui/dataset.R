dataset_tab <- tabItem(tabName="dataset",
  h1("Select an existing dataset or load a custom one"),
  fluidRow(
    box(
      title="Statistical parameters of the dataset",
      width=12,
      verbatimTextOutput("dataset_stats_info")
    ),
    box(
      title="Dataset table",
      width=12,
      dataTableOutput("dataset_table")
    )
  )
)

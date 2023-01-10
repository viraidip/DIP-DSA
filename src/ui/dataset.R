dataset_tab <- tabItem(tabName="dataset",
  h1("Overview about the selected dataset"),
  fluidRow(
    box(
      title="Statistical overview",
      width=12,
      "Different statistical parameters for the NGS count of the dataset",
      verbatimTextOutput("dataset_stats_info")
    ),
    box(
      title="Inspecting single entries",
      width=12,
      "Entries can be selected by clicking on them. Further information can",
      "then be seen in the",
      actionLink("link_to_single_datapoint_tab", "Inspect single datapoint"),
      " tab",
      dataTableOutput("dataset_table")
    )
  )
)

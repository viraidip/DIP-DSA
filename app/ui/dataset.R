dataset_tab <- tabItem(tabName="dataset",
  h1("Overview about the selected dataset"),
  fluidRow(
    splitLayout(
      cellWidths = c("50%", "50%"),
      box(
        title="Statistical overview",
        width=12,
        "Different statistical parameters for the NGS count of the dataset.",
        verbatimTextOutput("dataset_stats_info")
      ),
      box(
        title = "NGS count distribution",
        width = 12,
        plotlyOutput("ngs_distribution_plot")
      )
    ),
    box(
      title="Inspecting single entries",
      width=12,
      "Entries can be selected by clicking on them. Further information can.",
      "then be seen in the",
      actionLink("link_to_single_data_point_tab", "inspect single datapoint"),
      " tab",
      dataTableOutput("dataset_table")
    ),
    box(
      title="Venn diagram",
      width=12,
      "Displays a venn diagramm, if two data sets are selected. Shows how",
      "many DI candidates are in each set and how many of them can be found",
      "in both.",
      plotOutput("dataset_venn"),
      dataTableOutput("candidate_intersection_table")
    )
  )
)

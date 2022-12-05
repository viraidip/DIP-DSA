single_datapoint_tab <- tabItem(tabName="single_datapoint",
  h1("Inspect a single data point"),
  fluidRow(
    box(
      title="Info about the data point",
      width=12,
      "Go to ",
      actionLink("link_to_dataset_tab", "Data set overview tab"), 
      " to select a data point there at the bottom of the page.",
      br(),
      br(),
      tags$b("General info about the selected data point:"),
      verbatimTextOutput("single_datapoint_info"),
      br(),
      tags$b("Information about the packaging signal (if available):"),
      verbatimTextOutput("single_datapoint_packaging_signal_info")
    ),
    box(
      title="Nucleotides at start of deletion site",
      width=6,
      "Overview of the nucleotides directly at the start of the deletion",
      "site. The white area remains in the DI RNA sequence, the grey area is",
      "deleted.",
      plotOutput("single_datapoint_start_window")
    ),
    box(
      title="Nucleotides at end of deletion site",
      width=6,
      "Overview of the nucleotides directly at the end of the deletion site.",
      "The white area remains in the DI RNA sequence, the grey area is",
      "deleted.",
      plotOutput("single_datapoint_end_window")
    )
  )
)

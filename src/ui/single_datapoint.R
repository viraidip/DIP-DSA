single_datapoint_tab <- tabItem(tabName="single_datapoint",
  h1("Inspect a single data point"),
  fluidRow(
    box(
      title="Info about the data point",
      width=12,
      "Go to ",
      actionLink("link_to_dataset_tab", "load/select dataset tab"), 
      " to select a data point there at the bottom of the page.",
      br(),
      br(),
      verbatimTextOutput("single_datapoint_info")
    ),
    box(
      title="Start site",
      width=6,
      plotOutput("single_datapoint_start_window")
    ),
    box(
      title="End site",
      width=6,
      plotOutput("single_datapoint_end_window")
    )
  )
)

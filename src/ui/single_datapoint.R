single_datapoint_tab <- tabItem(tabName="single_datapoint",
  h1("Inspect a single datapoint"),
  fluidRow(
    box(
      title="Selected datapoint",
      width=12,
      "Go to ",
      actionLink("link_to_dataset_tab", "load/select dataset tab"), 
      " to select a data point there at the bottom of the page.",
      br(),
      br(),
      "The selected datapoint is row",
      uiOutput("single_datapoint")
    ),
    box(
      title="Start Site",
      width=6,
    ),
    box(
      title="End Site",
      width=6,
    )
  )
)

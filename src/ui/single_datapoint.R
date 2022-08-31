single_datapoint_tab <- tabItem(tabName="single_datapoint",
  h1("Inspect a single datapoint"),
  fluidRow(
    box(
      title="Selected datapoint",
      width=12,
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

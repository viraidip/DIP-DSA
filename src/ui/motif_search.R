regression_tab <- tabItem(tabName="motif_search",
  h1("Linear regression of segment length and NGS count"),
  fluidRow(
    box(
      title="Regression plot",
      width=12,
      selectInput(
        inputId="regression_segments",
        label="Select segments to include in linear regression",
        choices=SEGMENTS,
        selected=SEGMENTS,
        multiple=TRUE
      ),
      "The relative occurrence of each segment in the full dataset is",
      "calculated and plotted against the length of the segment. The expected",
      "values are the length of a single segment divided by the sum of all",
      "segment lengths.",
      plotOutput("regression_plot")
    )
  )
)

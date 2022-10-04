regression_tab <- tabItem(tabName="regression",
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
      plotOutput("regression_plot")
    )
  )
)

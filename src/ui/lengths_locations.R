lengths_locations_tab <- tabItem(tabName="lengths_locations",
  h1("Length and location of deletion sites"),
  fluidRow(
    box(
      title="Locations",
      width=12,
      selectInput(
        inputId="locations_segment",
        label="Select segment to show:",
        choices=SEGMENTS,
      ),
      radioButtons(
        inputId="locations_flattened",
        label="Show data flattened or unflattened (including NGS count)",
        choices=c("flattened", "unflattened"),
      ),
      plotOutput("locations_plot")
    ),
    box(
      title="Lengths",
      width=12,
      selectInput(
        inputId="lengths_segment",
        label="Select segment to show:",
        choices=SEGMENTS,
      ),
      radioButtons(
        inputId="lengths_flattened",
        label="Show data flattened or unflattened (including NGS count)",
        choices=c("flattened", "unflattened"),
      ),
      sliderInput(
        inputId="lengths_bins",
        label="Set size of bins for histogram",
        1, 100, 20
      ),
      plotOutput("lengths_plot")
    )
  )
)

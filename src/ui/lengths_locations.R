lengths_locations_tab <- tabItem(tabName="lengths_locations",
  h1("Length and location of deletion sites"),
  fluidRow(
    box(
      title="Lengths",
      width=12,
      selectInput(
        inputId="lengths_segment",
        label="Select segment to show:",
        choices=c("PB1"),
      ),
      radioButtons(
        inputId="lengths_flattened",
        label="Show data flattened or unflattened (including NGS count)",
        choices=c("flattened", "unflattened"),
      ),
      plotOutput("lengths_plot")
    )
  )
)

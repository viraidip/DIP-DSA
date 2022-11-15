lengths_locations_tab <- tabItem(tabName="lengths_locations",
  h1("Length and location of deletion sites"),
  fluidRow(
    box(
      title="Deletion site locations on full segment",
      width=12,
      radioButtons(
        inputId="locations_flattened",
        label="Show data flattened or unflattened (including NGS count)",
        choices=c("flattened", "unflattened"),
        inline=TRUE
      ),
      "Locations of all start and end positions of the deletion sites are",
      "plotted in reference to the full sequence.",
      plotlyOutput("locations_plot")
    ),
    box(
      title="Lengths of DI RNAs",
      width=12,
      radioButtons(
        inputId="lengths_flattened",
        label="Show data flattened or unflattened (including NGS count)",
        choices=c("flattened", "unflattened"),
        inline=TRUE
      ),
      sliderInput(
        inputId="lengths_bins",
        label="Set size of bins for histogram",
        1, 100, 20
      ),
      "The length of the single DI RNAs is plotted as a histogram, showing the",
      "number of occurrences for each length.",
      plotlyOutput("lengths_plot")
    )
  )
)

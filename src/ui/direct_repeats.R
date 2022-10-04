direct_repeats_tab <- tabItem(tabName="direct_repeats",
  h1("Direct Repeats around deletion sites"),
  fluidRow(
    box(
      title="direct repeats",
      width=12,
      radioButtons(
        inputId="direct_repeats_flattened",
        label="Flattened or unflattened data (including NGS count):",
        choices=c("flattened", "unflattened"),
        inline=TRUE
      ),
      plotlyOutput("direct_repeats_plot")
    )
  )
)

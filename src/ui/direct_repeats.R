direct_repeats_tab <- tabItem(tabName="direct_repeats",
  h1("Direct repeats around deletion sites"),
  fluidRow(
    box(
      title="Frequency of direct repeats",
      width=12,
      radioButtons(
        inputId="direct_repeats_flattened",
        label="Flattened or unflattened data (including NGS count):",
        choices=c("flattened", "unflattened"),
        inline=TRUE
      ),
      radioButtons(
        inputId="direct_repeats_correction",
        label="Perform adjustment/correction of data:",
        choices=c("Yes", "No"),
        inline=TRUE
      ),
      plotlyOutput("direct_repeats_plot")
    )
  )
)

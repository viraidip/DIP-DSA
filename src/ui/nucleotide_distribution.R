nucleotide_distribution_tab <- tabItem(tabName="nucleotide_distribution",
  h1("Nucleotide distribution at deletion sites"),
  fluidRow(
    box(
      title="Start site",
      width=6,
      radioButtons(
        inputId="nuc_dist_start_flattened",
        label="Show data flattened or unflattened (including NGS count)",
        choices=c("flattened", "unflattened"),
        inline=TRUE
      ),
      plotOutput("nuc_dist_start")
    ),
    box(
      title="End site",
      width=6,
      radioButtons(
        inputId="nuc_dist_end_flattened",
        label="Show data flattened or unflattened (including NGS count)",
        choices=c("flattened", "unflattened"),
        inline=TRUE
      ),
      plotOutput("nuc_dist_end")
    )
  )
)

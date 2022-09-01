nucleotide_distribution_tab <- tabItem(tabName="nucleotide_distribution",
  h1("Nucleotide distribution at deletion sites"),
  fluidRow(
    box(
      title="Start site",
      width=6,
      plotOutput("nuc_dist_start")
    ),
    box(
      title="End site",
      width=6,
      plotOutput("nuc_dist_end")
    )
  )
)

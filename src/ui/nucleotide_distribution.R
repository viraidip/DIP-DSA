nucleotide_distribution_tab <- tabItem(tabName="nucleotide_distribution",
  h1("Nucleotide distribution at deletion sites"),
  fluidRow(
    box(
      title="Set parameters",
      width=12,
      radioButtons(
        inputId="nuc_dist_flattened",
        label="Show data flattened or unflattened (including NGS count)",
        choices=c("flattened", "unflattened"),
        inline=TRUE
      )
    ),
    box(
      title="Start site",
      width=6,
      plotOutput("nuc_dist_start_A"),
      plotOutput("nuc_dist_start_C"),
      plotOutput("nuc_dist_start_G"),
      plotOutput("nuc_dist_start_U")
    ),
    box(
      title="End site",
      width=6,
      plotOutput("nuc_dist_end_A"),
      plotOutput("nuc_dist_end_C"),
      plotOutput("nuc_dist_end_G"),
      plotOutput("nuc_dist_end_U")
    )
  )
)

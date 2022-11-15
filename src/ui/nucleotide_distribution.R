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
      ),
      "The relative occurrence for each of the four nucleotides is calculated",
      "(observed) and compared to a random sampling (expected)."
    ),
    box(
      title="Start of deletion site",
      width=6,
      "Distribution of the four nucleotides at the start of the deletion site.",
      plotlyOutput("nuc_dist_start_A"),
      plotlyOutput("nuc_dist_start_C"),
      plotlyOutput("nuc_dist_start_G"),
      plotlyOutput("nuc_dist_start_U")
    ),
    box(
      title="End of deletion site",
      width=6,
      "Distribution of the four nucleotides at the end of the deletion site.",
      plotlyOutput("nuc_dist_end_A"),
      plotlyOutput("nuc_dist_end_C"),
      plotlyOutput("nuc_dist_end_G"),
      plotlyOutput("nuc_dist_end_U")
    )
  )
)

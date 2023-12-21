multiple_datasets_tab <- tabItem(tabName="multiple_datasets",
  h1("Analyse multiple datasets"),
  fluidRow(
    box(
      title="Select datasets",
      width=12,
      selectInput(
        inputId="multiple_datasets",
        label="Select datasets to compare:",
        choices=list.files(DATASETSPATH,
          pattern="csv$",
          full.names=FALSE,
          recursive=TRUE
        ),
        multiple=TRUE,
        selected=c("A_California_07_2009/Alnaji2019_Cal07.csv", "A_NewCaledonia_1999/Alnaji2019_NC.csv")
      ),
      sliderInput(
        inputId="multiple_RCS",
        label="Set RCS:",
        1,
        100,
        2,
        step=1
      ),
      radioButtons(
        inputId="multiple_flattened",
        label="Show data flattened or unflattened (including NGS count):",
        choices=c("flattened", "unflattened"),
        inline=TRUE
      ),
      selectInput(
        inputId="multiple_selected_segment",
        label="Select segment",
        choices=c(SEGMENTS, "ALL")
      ),
    ),
    box(
      title="NGS counts",
      width=12,
      "Different statistical parameters for the NGS count of the datasets",
      plotlyOutput("multiple_ngs_distribution_plot")
    ),
    box(
      title="Deletion shift",
      width=12,
      "Distribution of the deletion shift for the selected datasets",
      plotlyOutput("multiple_deletion_shift_plot")
    ),
    box(
      title="Segment distribution",
      width=12,
      "Distribution on the eight RNA segments of the selected datasets",
      plotlyOutput("multiple_segment_distribution_plot")
    ),
    box(
      title="DVG length",
      width=12,
      sliderInput(
        inputId="multiple_lengths_bins",
        label="Set size of bins for histogram:",
        1,
        1000,
        20,
        step=1
      ),
      "Lengths of the DVGs",
      plotlyOutput("multiple_deletion_length_plot")
    ),
    box(
      title="Nucleotide enrichment (Start position)",
      width=6,
      "Distribution of nucleotide enrichment at start of the deletion site.",
      radioButtons(
        inputId="multiple_enrichment_nucleotide_start",
        label="Select nucleotide:",
        choices=c("A", "C", "G", "U"),
        inline=TRUE
      ),
      plotlyOutput("multiple_nucleotide_enrichment_start_plot")
    ),
    box(
      title="Nucleotide enrichment (End position)",
      width=6,
      "Distribution of nucleotide enrichment at end of the deletion site.",
      radioButtons(
        inputId="multiple_enrichment_nucleotide_end",
        label="Select nucleotide:",
        choices=c("A", "C", "G", "U"),
        inline=TRUE
      ),
      plotlyOutput("multiple_nucleotide_enrichment_end_plot")
    ),
    box(
      title="Direct repeats",
      width=12,
      "Distribution of the direct repeat lengths.",
      plotlyOutput("multiple_direct_repeats_plot")
    )
  )
)



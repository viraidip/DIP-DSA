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
        selected=c("A_California_07_2009/Alnaji2019.csv", "A_PuertoRico_8_1934/alnaji2021_extra.csv")
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
      title="Segment distribution",
      width=12,
      "Distribution on the eight RNA segments of the selected datasets",
      plotlyOutput("multiple_segment_distribution_plot")
    ),

    box(
      title="Deletion shift",
      width=12,
      "Distribution of the deletion shift for the selected datasets",
      plotlyOutput("multiple_deletion_shift_plot")
    ),

    box(
      title="DVG length",
      width=12,
      "Lengths of the DVGs",
      plotlyOutput("multiple_deletion_length_plot")
    ),

    box(
      title="Nucleotide enrichment",
      width=12,
      "Distribution of the deletion shift for the selected datasets",
      plotlyOutput("multiple_nucleotide_enrichment_plot")
    ),

    box(
      title="Direct repeats",
      width=12,
      "Distribution of the deletion shift for the selected datasets",
      plotlyOutput("multiple_direct_repeats_plot")
    ),





  )
)



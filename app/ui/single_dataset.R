single_dataset_tab <- tabItem(tabName="single_dataset",
  h1("Analyse single dataset"),
  fluidRow(
    box(
      width=12,
      title="Select data",
      selectInput(
        inputId="single_strain",
        label="Select a strain:",
        choices=gsub(
          "_",
          "/",
          list.dirs(DATASETSPATH, full.names=FALSE, recursive=FALSE)
        )
      ),
      selectInput(
        inputId="single_dataset",
        label="Select an existing dataset:",
        choices="Alnaji2019_Cal07"
      ),
      sliderInput(
        inputId="single_RSC",
        label="Set the RSC (read support cutoff):",
        1,
        100,
        2,
        step=1
      ),
      radioButtons(
        inputId="single_flattened",
        label="Show data flattened or unflattened (including NGS count):",
        choices=c("flattened", "unflattened"),
        inline=TRUE
      ),
      selectInput(
        inputId="single_selected_segment",
        label="Select segment:",
        choices=c(SEGMENTS)
      ),
      actionButton(
        inputId="single_submit",
        label="Generate plots"
      )
    ),
    box(
      title="Distribution of NGS count",
      width=6,
      "The NGS counts of the DVGs given in the selected dataset.",
      plotlyOutput("ngs_distribution_plot")
    ),
    box(
      title="Deletion shift",
      width=6,
      "Shift of the reading frame introduced by deletion site.",
      plotlyOutput("deletion_shift_plot"),
    ),
    box(
      title="Segment distribution",
      width=6,
      "Distribution of the DVGs over the eight segments.",
      plotlyOutput("segment_distribution_plot"),
    ),
    box(  
      title="DVG lengths",
      width=6,
      sliderInput(
        inputId="single_lengths_bins",
        label="Set size of bins for histogram:",
        1,
        400,
        20,
        step=1
      ),
      "The length of the single DVGs is plotted as a histogram, showing",
      "the number of occurrences for each length.",
      plotlyOutput("lengths_plot"),
    ),
    box(
      title="Deletion locations",
      width=6,
      "Locations and nucleotides of all start and end positions of the ",
      "deletion sites are plotted in reference to the full sequence.",
      plotlyOutput("locations_plot"),
    ),
    box(
      title="3' and 5' sequence end comparision",
      width=6,
      "Comparision of the lengths of the 3' and 5' ends. If data about the",
      "packaging signal is available, incorporation signal is included in",
      "blue and bundling signal is included in red",
      plotlyOutput("end_3_5_plot")
    ),
    box(
      title="Deletion connections",
      width=6,
      "The plot displays the connection between the start and end positons",
      "of the DVGs.",
      plotlyOutput("start_end_connection_plot"),
    ),

    box(
      title="Frequency of direct repeats",
      width=6,
      "The length of the overlapping sequence of the start and end of the",
      "deletion site is calculated and plotted in a bar plot. The results are",
      "compared against a sampling apporach.",
      plotlyOutput("direct_repeats_plot")
    ),
    box(
      title="Nucleotide enrichment at start of deletion site",
      width=6,
      radioButtons(
        inputId="enrichment_nucleotide_start",
        label="Select nucleotide:",
        choices=c("A", "C", "G", "U"),
        inline=TRUE
      ),
      plotlyOutput("nucleotide_enrichment_start_plot"),
    ),
    box(
      title="Nucleotide enrichment at end of deletion site",
      width=6,
      radioButtons(
        inputId="enrichment_nucleotide_end",
        label="Select nucleotide:",
        choices=c("A", "C", "G", "U"),
        inline=TRUE
      ),
      plotlyOutput("nucleotide_enrichment_end_plot")
    )
  )
)

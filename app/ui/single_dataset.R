single_dataset_tab <- tabItem(tabName="single_dataset",
  h1("Analyse single dataset"),
  fluidRow(
    box(
      width=12,
      title="Select data",
      "More info about the predefined datasets can be found in the",
      actionLink("link_to_about_tab", "about"),
      "tab.",
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
        choices="Alnaji2019"
      ),
      sliderInput(
        inputId="single_RCS",
        label="Set RCS:",
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
        label="Select segment",
        choices=c(SEGMENTS, "ALL")
      ),

    ),
    box(
      title="Statistical overview",
      width=6,
      "Different statistical parameters for the NGS count of the dataset.",
      plotlyOutput("ngs_distribution_plot")
    ),

    # deletion shift
    box(
      title="Deletion shift",
      width=6,
      "Shift of the reading frame introduced by deletion site.",
      plotlyOutput("deletion_shift_plot"),
    ),

    # segment distribution
    box(
      title="Segment distribution",
      width=12,
      "Distribution of the DVGs over the eight segments.",
      plotlyOutput("segment_distribution_plot"),
    ),

    # deletion length and location
    box(  
      title="DVG lengths",
      width=12,
      sliderInput(
        inputId="single_lengths_bins",
        label="Set size of bins for histogram:",
        1,
        100,
        20
      ),
      "The length of the single DI RNAs is plotted as a histogram, showing",
      "the number of occurrences for each length.",
      plotlyOutput("lengths_plot"),

      "Locations of all start and end positions of the deletion sites are",
      "plotted in reference to the full sequence.",
      plotlyOutput("locations_plot"),
      "Descriptive text",
      plotlyOutput("start_end_connection_plot"),
    
    
      "Comparision of the lengths of the 3' and 5' ends. If data about the",
      "packaging signal is available, incorporation signal is included in",
      "blue and bundling signal is included in red",
      plotlyOutput("end_3_5_plot")
    ),

    # nucleotide enrichment
    box(
      title="Nucleotide enrichment at start of deletion site",
      width=6,
      "Distribution of the four nucleotides at the start of the deletion site.",
      plotlyOutput("nuc_dist_start_A"),
      plotlyOutput("nuc_dist_start_C"),
      plotlyOutput("nuc_dist_start_G"),
      plotlyOutput("nuc_dist_start_U")
    ),
    box(
      title="Nucleotide enrichment at end of deletion site",
      width=6,
      "Distribution of the four nucleotides at the end of the deletion site.",
      plotlyOutput("nuc_dist_end_A"),
      plotlyOutput("nuc_dist_end_C"),
      plotlyOutput("nuc_dist_end_G"),
      plotlyOutput("nuc_dist_end_U")
    ),

    # direct repeats
    box(
      title="Frequency of direct repeats",
      width=12,
      "The length of the overlapping sequence of the start and end of the",
      "deletion site is calculated and plotted in a bar plot. The results are",
      "compared against a sampling apporach.",
      plotlyOutput("direct_repeats_plot")
    ),

  )
)

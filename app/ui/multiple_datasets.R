multiple_datasets_tab <- tabItem(tabName="multiple_datasets",
  h1("Analyse multiple datasets"),
  fluidRow(
    box(
      title="Select dataset and define parameters for analysis",
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
        selected=c("A_California_07_2009/Alnaji2019_Cal07.csv",
          "A_NewCaledonia_20-JY2_1999/Alnaji2019_NC.csv"
        )
      ),
      sliderInput(
        inputId="multiple_RSC",
        label="Set the RSC (read support cutoff):",
        1,
        100,
        2,
        step=1
      ),
      "The RSC is dataset-specific and is usually set to a value between 5",
      "and 30. As reference: In our meta-analysis an RSC of 15 was used.",
      br(),
      br(),
      radioButtons(
        inputId="multiple_flattened",
        label="Show data flattened or unflattened (including NGS count):",
        choices=c("flattened", "unflattened"),
        inline=TRUE
      ),
      "Flattened data do not take the NGS count into account. In unflattened",
      "data, each individual DelVG is weighted by the NGS count.",
      br(),
      br(),
      selectInput(
        inputId="multiple_selected_segment",
        label="Select segment:",
        choices=c(SEGMENTS)
      ),
      actionButton(
        inputId="multiple_submit",
        label="Generate plots"
      )
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
      "Distribution of the reading frame shift for the selected datasets. A",
      "chi-squared test is performed to compare the distribution against a",
      "random shifts.",
      plotlyOutput("multiple_deletion_shift_plot")
    ),
    box(
      title="Segment distribution",
      width=12,
      "Distribution of the DelVGs over the eight segments. It is tested by a",
      "chi-squared test if the distribution is similar to a distribution that",
      "would be expected if the DelVGs occur solely dependent on the RNA",
      "sequence length.",
      plotlyOutput("multiple_segment_distribution_plot")
    ),
    box(
      title="DelVG length",
      width=12,
      sliderInput(
        inputId="multiple_lengths_bins",
        label="Set size of bins for histogram:",
        1,
        1000,
        20,
        step=1
      ),
      "Length distribution of the DelVGs. Datasets can be (un-)selected by",
      "clicking on them in the legend.",
      plotlyOutput("multiple_deletion_length_plot")
    ),
    box(
      title="Nucleotide enrichment (Start position)",
      width=6,
      "Distribution of nucleotide enrichment at start of the deletion site.",
      "For each position the difference to randomly sampled data is estimated",
      "using a one-way ANOVA.",
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
      "For each position the difference to randomly sampled data is estimated",
      "using a one-way ANOVA.",
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
      "Distribution of the direct repeat lengths. The distribution is",
      "compared against data from a random sampling apporach using a chi",
      "squared test.",
      plotlyOutput("multiple_direct_repeats_plot")
    )
  )
)



dataset_intersection_tab <- tabItem(tabName="dataset_intersection",
  h1("Overview about all dataset"),
  fluidRow(
    box(
      title="Intersecting candidates",
      width=12,
      selectInput(
        inputId="selected_datasets",
        label="Select datasets to compare:",
        choices=list.files(DATASETSPATH,
          pattern="csv$",
          full.names=FALSE,
          recursive=TRUE
        ),
        multiple=TRUE,
        selected=c("A_California_07_2009/Alnaji2019_Cal07.csv",
          "A_NewCaledonia_1999/Alnaji2019_NC.csv"
        )
      ),
      sliderInput(
        inputId="RSC_intersection",
        label="Set the RSC (read support cutoff)",
        1,
        100,
        2,
        step=1
      ),
      actionButton(
        inputId="intersection_submit",
        label="Generate plots"
      )
    ),
    box(
      title="DelVG candidates above defined threshold",
      width=12,
      sliderInput(
        inputId="min_occurrences",
        label="Set the number of minimum occurrences of a DelVG candidate:",
        1,
        length(list.files(DATASETSPATH,
          pattern="csv$",
          full.names=FALSE,
          recursive=TRUE
          )
        ),
        2,
        step=1
      ),
      dataTableOutput("intersecting_candidates_table")
    ),
    box(
      title="Intersecting DelVG candidates",
      width=12,
      plotlyOutput("overlap_matrix_plot")
    ),
    box(
      title="Upset plot",
      width=12,
      plotOutput("upset_plot")
    ),
    box(
      title="Identified candidates with NGS counts",
      width=12,
      plotlyOutput("intersecting_candidates_NGS_plot")
    )
  )
)

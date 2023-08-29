all_datasets_tab <- tabItem(tabName="all_datasets",
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
        multiple=TRUE
      ),
      sliderInput(
        inputId="min_occurrences",
        label="Set the number of minimum occurrences of a DI candidate:",
        1,
        length(list.files(DATASETSPATH,
          pattern="csv$",
          full.names=FALSE,
          recursive=TRUE
          )
        ),
        2,
        step=1
      )
    ),
    box(
      title="table",
      width=12,
      dataTableOutput("intersecting_candidates_table")
    ),
    box(
      title="Intersecting DI candidates",
      width=12,
      plotlyOutput("overlap_matrix_plot")
    ),
    box(
      title="Upset plot",
      width=12,
      plotOutput("upset_plot")
    )
  )
)

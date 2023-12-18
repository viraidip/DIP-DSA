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
        selected=c("A_California_07_2009/Alnaji2019.csv", "A_PuertoRico_8_1934/alnaji2021_extra.csv")
      ),
      sliderInput(
        inputId="RCS_intersection",
        label="Set RCS:",
        1,
        100,
        2,
        step=1
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

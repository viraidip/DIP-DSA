dataset_intersection_tab <- tabItem(tabName="dataset_intersection",
  h1("Identification of promising DelVGs by intersecting multiple datasets"),
  fluidRow(
    box(
      title="Select datasets for intersection analysis",
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
      actionButton(
        inputId="intersection_submit",
        label="Generate plots"
      )
    ),
    box(
      title="Intersecting DelVG candidates per dataset pair",
      width=12,
      plotlyOutput("overlap_matrix_plot")
    ),
    box(
      title="Number of candidates per segment that occur in multiple datasets",
      width=12,
      plotlyOutput("barplot_candidates_plot")
    ),
    box(
      title="Identified candidates by using mean and sum scoring",
      width=12,
      selectInput(
        inputId="intersection_selected_segment",
        label="Select segment:",
        choices=c(SEGMENTS)
      ),
      sliderInput(
        inputId="intersection_thresh",
        label="Using the highest n ranked candidates:",
        1,
        100,
        10,
        step=1
      ),
      plotlyOutput("highest_n_ranked_plot")
    )
  )
)

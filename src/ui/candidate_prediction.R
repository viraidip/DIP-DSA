candidate_prediction_tab <- tabItem(tabName="candidate_prediction",
  h1("Predicting the competitveness of a DI RNA candidate"),
  fluidRow(
    box(
      title="Defining DI RNA candidate to predict",
      width=12,
      numericInput(
        inputId="prediction_start",
        label="Type the starting position:",
        value=100,
        min=0,
        max=1000
      ),
      numericInput(
        inputId="prediction_end",
        label="Type the ending position:",
        value=400,
        min=0,
        max=1000
      ),
      selectInput(
        inputId="prediction_segment",
        label="Select a segment:",
        choices=SEGMENTS
      ),
      selectInput(
        inputId="prediction_strain",
        label="Select a strain:",
        choices=gsub(
          "_",
          "/",
          list.dirs(DATASETSPATH, full.names=FALSE, recursive=FALSE)
        )
      ),
      selectInput(
        inputId="classifier",
        label="Select a classifier to use:",
        choices=gsub(
          ".pkl",
          "",
          list.files(
            file.path(DATAPATH, "classifiers"),
            full.names=FALSE,
            recursive=FALSE
          )
        )
      ),
      actionButton(
        inputId="run_prediction",
        label="Run prediction"
      )
    ),
    box(
      title="Prediction result",
      width=12,
      "Displaying the result of the prediction. If no result is shown press",
      "'Run prediction'. The calculations can take a few minutes.",
      verbatimTextOutput("prediction_results")
    )
  )
)

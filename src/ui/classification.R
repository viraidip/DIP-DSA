classification_tab <- tabItem(tabName="classification",
  h1("Predicting the competitveness of a DI RNA candidate"),
  fluidRow(
    box(
      title="Defining DI RNA candidate to classify",
      width=12,
      numericInput(
        inputId="clf_start",
        label="Define start position:",
        value=100
      ),
      numericInput(
        inputId="clf_end",
        label="Define end position:",
        value=800
      ),
      selectInput(
        inputId="clf_segment",
        label="Select segment:",
        choices=SEGMENTS
      ),
      selectInput(
        inputId="clf_strain",
        label="Select a strain:",
        choices=gsub(
          "_",
          "/",
          list.dirs(DATASETSPATH, full.names=FALSE, recursive=FALSE)
        )
      ),
      selectInput(
        inputId="classifier",
        label="Select classifier",
        choices=c("clf A", "clf B")
      ),
      actionButton(
        inputId="run_clf",
        label="Run classification"
      )
    ),
    box(
      title="Classification result",
      width=12,
      "text"
    )
  )
)

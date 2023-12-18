multiple_datasets_tab <- tabItem(tabName="multiple_datasets",
  h1("Analyse multiple datasets"),
  fluidRow(
    box(
      title="Select datasets",
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
        inputId="mRCS_multiple",
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

  )
)

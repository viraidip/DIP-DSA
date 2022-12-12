np_density_tab <- tabItem(tabName="np_density",
  h1("Nucleoprotein density"),
  fluidRow(
    box(
      title="Input",
      width=12,
      "NP areas need to be defined in a single string. Each area is split by",
      "a comma and the start and end of each area is seperated by a hyphen",
      "(e.g. '10-40,60-85').",
      textInput(
        inputId="np_areas",
        label="Write the areas of high NP density"
      )
    ),
    box(
      title="Input",
      width=12,
      plotlyOutput("np_plot")
    ),
    box(
      title="Input",
      width=12,
      plotlyOutput("np_bar_plot")
    )
  )
)

np_density_tab <- tabItem(tabName="np_density",
  h1("Nucleoprotein density"),
  fluidRow(
    box(
      title="Input NP high areas",
      width=12,
      "NP areas need to be defined as a string. The start and end point of an",
      "area is split by a hyphen. To define multiple NP high areas seperate","
      them by a comma (e.g. '10-40,60-85').",
      textInput(
        inputId="np_areas",
        label="Write areas high NP density",
        placeholder="45-178,406-524,1676-1894"
      )
    ),
    box(
      title="Mapping of NP high areas to deletion sites",
      width=12,
      plotlyOutput("np_plot")
    ),
    box(
      title="Comparision of distribution",
      width=12,
      plotlyOutput("np_bar_plot")
    )
  )
)

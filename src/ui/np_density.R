np_density_tab <- tabItem(tabName="np_density",
  h1("Nucleoprotein density"),
  fluidRow(
    box(
      title="Input areas",
      width=12,
      "NP areas need to be defined as a string. The start and end point of an",
      "area is split by a hyphen. To define multiple NP high areas seperate","
      them by a comma (e.g. '10-40,60-85').",
      textInput(
        inputId="np_areas",
        label="Write areas of high NP density:",
        placeholder="45-178,406-524,1676-1894"
      )
    ),
    box(
      title="Mapping of high NP areas to deletion sites",
      width=12,
      "The plot shows the high NP areas in correspondence to the start and",
      "end positions of the deletion sites of the data set",
      plotlyOutput("np_plot")
    ),
    box(
      title="Comparision of distribution",
      width=12,
      "The number of deletion sites in high NP areas is divided by the",
      "overall count of deletion sites. The ratio of the observed samples is",
      "compared to randomly sampled data.",
      plotlyOutput("np_bar_plot")
    )
  )
)

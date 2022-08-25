
ui <- dashboardPage(
  dashboardHeader(title="DIP DSA"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Select/Load Dataset", tabName="dataset", icon=icon("database")),
      hr(),
      menuItem("Lengths and Locations", tabName="lengths_locations", icon=icon("ruler-horizontal")),
      menuItem("Nucleotide Distribution", tabName="nucleotide_distribution", icon=icon("magnifying-glass-chart")),
      menuItem("Direct Repeats", tabName="direct_repeats", icon=icon("repeat")),
      menuItem("Linear Regression", tabName="regression", icon=icon("chart-line")),
      hr(),
      menuItem("About", tabName="about", icon=icon("info"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName="dataset",
        h1("Select/Load dataset"),
        fluidRow(
        )
      ),
      tabItem(tabName="lengths_locations",
        h1("Length and location of deletion sites"),
        fluidRow(
        )
      ),
      tabItem(tabName="nucleotide_distribution",
        h1("Nucleotide distribution at deletion sites"),
        fluidRow(
        )
      ),
      tabItem(tabName="direct_repeats",
        h1("Direct Repeats around deletion sites"),
        fluidRow(
        )
      ),
      tabItem(tabName="regression",
        h1("Linear regression of segment length and NGS count"),
        fluidRow(
        )
      ),
      tabItem(tabName="about",
        h1("About this application"),
        fluidRow(
        )
      )
    )
  )
)


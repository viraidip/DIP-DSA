library(shiny)
library(shinydashboard)
library(ggplot2)


### UI ###
# Loads the sources for the UI of each tab.
# Each tab is saved in an individual file.
source('ui/dataset.R', local=TRUE)
source('ui/lengths_locations.R', local=TRUE)
source('ui/nucleotide_distribution.R', local=TRUE)
source('ui/direct_repeats.R', local=TRUE)
source('ui/regression.R', local=TRUE)
source('ui/about.R', local=TRUE)

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
      dataset_tab,
      lengths_locations_tab,
      nucleotide_distribution_tab,
      direct_repeats_tab,
      regression_tab,
      about_tab
    )
  )
)


### SERVER ###
# Load the sources for the server logic.
# Each tab has an own file for its server functions.
source("server/dataset.R", local=TRUE)
source("server/lengths_locations.R", local=TRUE)
#source("server/nucleotide_distribution.R", local=TRUE)
#source("server/direct_repeats.R", local=TRUE)
#source("server/regression.R", local=TRUE)
#source("server/about.R", local=TRUE)

server <- function(input, output) {
  ### load/select dataset ###
  load_dataset <- reactive({
    path <- file.path("..", "data", paste(input$dataset_name, ".csv", sep=""))
    read.csv(path)
  })

  output$dataset_table <- renderTable(load_dataset())

  ### lenghts and locations

  output$lengths_plot <- renderPlot({
    create_lengths_plot(load_dataset(), input$lengths_segment, input$lengths_flattened)
  })
  
}


### APP ###
# Main function that builds the application
shinyApp(ui=ui, server=server)

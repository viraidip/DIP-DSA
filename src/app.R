library(shiny)
library(shinydashboard)
library(ggplot2)
library(DT)
library(shinyWidgets)

source("utils.R")

##########
### UI ###
##########
# Loads the sources for the UI of each tab.
# Each tab is saved in an individual file.
source('ui/dataset.R', local=TRUE)
source('ui/lengths_locations.R', local=TRUE)
source('ui/nucleotide_distribution.R', local=TRUE)
source('ui/direct_repeats.R', local=TRUE)
source('ui/regression.R', local=TRUE)
source('ui/single_datapoint.R', local=TRUE)
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
      menuItem("Inspect single datapoint", tabName="single_datapoint", icon=icon("magnifying-glass")),
      hr(),
      menuItem("About", tabName="about", icon=icon("info"))
    ),
    selectInput(
      inputId="selected_segment",
      label="Select segment",
      choices=SEGMENTS
    )
  ),
  dashboardBody(
    tabItems(
      dataset_tab,
      lengths_locations_tab,
      nucleotide_distribution_tab,
      direct_repeats_tab,
      regression_tab,
      single_datapoint_tab,
      about_tab
    )
  )
)

##############
### SERVER ###
##############
# Load the sources for the server logic.
# Each tab has an own file for its server functions.
source("server/dataset.R", local=TRUE)
source("server/lengths_locations.R", local=TRUE)
source("server/nucleotide_distribution.R", local=TRUE)
#source("server/direct_repeats.R", local=TRUE)
#source("server/regression.R", local=TRUE)
#source("server/single_datapoint.R", local=TRUE)
#source("server/about.R", local=TRUE)

server <- function(input, output) {
  ### load/select dataset ###
  load_dataset <- reactive({
    path <- file.path("..", "data", "datasets", paste(input$strain, ".csv", sep=""))
    read.csv(path, na.strings=c("NaN"))
  })

  output$dataset_table <- renderDataTable(
    datatable(load_dataset(),
      options = list(
        pageLength=100
      ),
      selection="single"
    )
  )


  ### lenghts and locations ###
  output$locations_plot <- renderPlot({
    create_locations_plot(load_dataset(),
      input$selected_segment,
      input$locations_flattened
    )
  })

  output$lengths_plot <- renderPlot({
    create_lengths_plot(load_dataset(),
      input$selected_segment,
      input$strain,
      input$lengths_flattened,
      input$lengths_bins
    )
  })


  ### nucleotide distribution ###
  output$nuc_dist_start <- renderPlot({
    create_nuc_dist_plot(load_dataset(),
      input$strain,
      input$selected_segment,
      "Start",
      input$nuc_dist_start_flattened
    )
  })
  
  output$nuc_dist_end <- renderPlot({
    create_nuc_dist_plot(load_dataset(),
      input$strain,
      input$selected_segment,
      "End",
      input$nuc_dist_end_flattened
    )
  })


  ### single datapoint ###
  output$single_datapoint <- renderUI({
    input$dataset_table_rows_selected
  })

}

###########
### APP ###
###########
# Main function that builds the application
shinyApp(ui=ui, server=server)

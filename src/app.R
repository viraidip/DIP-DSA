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
  skin="red",
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

server <- function(input, output, session) {
  ### load/select dataset ###
  load_dataset <- reactive({
    path <- file.path(DATASETSPATH, paste(input$strain, ".csv", sep=""))
    read.csv(path, na.strings=c("NaN"))
  })

  observeEvent(input$dataset_submit, {
    # check if all fields are filled
    req(
      input$upload_strain, input$upload_dataset_file,
      input$upload_PB2_file, input$upload_PB1_file,
      input$upload_PA_file, input$upload_HA_file,
      input$upload_NP_file, input$upload_NA_file,
      input$upload_M_file, input$upload_NS_file
    )

    # move submitted files to right folder
    from_list <- list(
      input$upload_dataset_file$datapath,
      input$upload_PB2_file$datapath, input$upload_PB1_file$datapath,
      input$upload_PA_file$datapath, input$upload_HA_file$datapath,
      input$upload_NP_file$datapath, input$upload_NA_file$datapath,
      input$upload_M_file$datapath,input$upload_NS_file$datapath
    )
    to_list <- list(file.path(DATASETSPATH, paste(input$upload_strain, ".csv", sep="")))
    fasta_path <- file.path(FASTAPATH, input$upload_strain)
    dir.create(fasta_path)
    for (s in SEGMENTS) {
      to_list <- append(to_list, file.path(fasta_path, paste(s, ".fasta", sep="")))
    }
    move_files(from_list, to_list)

    # select the new submitted dataset
    all_datasets <- tools::file_path_sans_ext(list.files(DATASETSPATH))
    updateSelectInput(
      session, inputId="strain", label="strain",choices=all_datasets,
      selected=input$upload_strain
      )
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

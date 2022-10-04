library(shiny)
library(shinydashboard)
library(ggplot2)
library(DT)
library(shinyWidgets)
library(plotly)

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

ui <- bootstrapPage(
  dashboardPage(
    skin="red",
    dashboardHeader(title="DIP DSA"),
    dashboardSidebar(
      sidebarMenu(
        id="sidebarmenu",
        menuItem("Select/Load Dataset", tabName="dataset", icon=icon("database")),
        menuItem("Inspect single datapoint", tabName="single_datapoint", icon=icon("magnifying-glass")),
        hr(),
        selectInput(
          inputId="selected_segment",
          label="Select segment",
          choices=SEGMENTS
        ),
        menuItem("Lengths and Locations", tabName="lengths_locations", icon=icon("ruler-horizontal")),
        menuItem("Nucleotide Distribution", tabName="nucleotide_distribution", icon=icon("magnifying-glass-chart")),
        menuItem("Direct Repeats", tabName="direct_repeats", icon=icon("repeat")),
        hr(),
        menuItem("Linear Regression", tabName="regression", icon=icon("chart-line")),
        hr(),
        menuItem("About", tabName="about", icon=icon("info"))
      )
    ),
    dashboardBody(
      tabItems(
        dataset_tab,
        single_datapoint_tab,
        lengths_locations_tab,
        nucleotide_distribution_tab,
        direct_repeats_tab,
        regression_tab,
        about_tab
      )
    )
  )
)

##############
### SERVER ###
##############
# Load the sources for the server logic.
# Each tab has an own file for its server functions.
source("server/dataset.R", local=TRUE)
source("server/single_datapoint.R", local=TRUE)
source("server/lengths_locations.R", local=TRUE)
source("server/nucleotide_distribution.R", local=TRUE)
source("server/direct_repeats.R", local=TRUE)
source("server/regression.R", local=TRUE)
#source("server/about.R", local=TRUE)

server <- function(input, output, session) {
### load/select dataset ###
  load_dataset <- reactive({
    path <- file.path(DATASETSPATH, paste(input$strain, ".csv", sep=""))
    col_names <- c("Segment", "Start", "End", "NGS_read_count")
    col_classes <- c("character", "integer", "integer", "integer")
    read.csv(path, na.strings=c("NaN"), col.names=col_names, colClasses=col_classes)
  })

  observeEvent(input$link_to_about_tab, {
    updateTabItems(session, "sidebarmenu", "about")
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
      input$upload_M_file$datapath, input$upload_NS_file$datapath
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


### single datapoint ###
  observeEvent(input$link_to_dataset_tab, {
    updateTabItems(session, "sidebarmenu", "dataset")
  })

  output$single_datapoint_info <- renderText({
    create_single_datapoint_info(load_dataset(),
      input$dataset_table_rows_selected,
      input$strain
    )
  })

  output$single_datapoint_start_window <- renderPlot({
    plot_deletion_site_window(load_dataset(),
      input$dataset_table_rows_selected,
      input$strain,
      "Start"
    )
  })

  output$single_datapoint_end_window <- renderPlot({
    plot_deletion_site_window(load_dataset(),
      input$dataset_table_rows_selected,
      input$strain,
      "End"
    )
  })


### lenghts and locations ###
  output$locations_plot <- renderPlotly({
    create_locations_plot(load_dataset(),
      input$selected_segment,
      input$locations_flattened
    )
  })

  output$lengths_plot <- renderPlotly({
    create_lengths_plot(load_dataset(),
      input$selected_segment,
      input$strain,
      input$lengths_flattened,
      input$lengths_bins
    )
  })


### nucleotide distribution ###
  observeEvent(input$strain, {
    create_nuc_dist_data(load_dataset(),
      input$strain,
      input$selected_segment,
      input$nuc_dist_flattened
    )
    update_plots()
  })
  observeEvent(input$selected_segment, {
    create_nuc_dist_data(load_dataset(),
      input$strain,
      input$selected_segment,
      input$nuc_dist_flattened
    )
    update_plots()
  })
  observeEvent(input$nuc_dist_flattened, {
    create_nuc_dist_data(load_dataset(),
      input$strain,
      input$selected_segment,
      input$nuc_dist_flattened
    )
    update_plots()
  })

  update_plots <- function() {
    output$nuc_dist_start_A <- renderPlotly({
      create_nuc_dist_plot("Start", "A")
    })
    output$nuc_dist_start_C <- renderPlotly({
      create_nuc_dist_plot("Start", "C")
    })
    output$nuc_dist_start_G <- renderPlotly({
      create_nuc_dist_plot("Start", "G")
    })
    output$nuc_dist_start_U <- renderPlotly({
      create_nuc_dist_plot("Start", "U")
    })
  
    output$nuc_dist_end_A <- renderPlotly({
      create_nuc_dist_plot("End", "A")
    })
    output$nuc_dist_end_C <- renderPlotly({
      create_nuc_dist_plot("End", "C")
    })
    output$nuc_dist_end_G <- renderPlotly({
      create_nuc_dist_plot("End", "G")
    })
    output$nuc_dist_end_U <- renderPlotly({
      create_nuc_dist_plot("End", "U")
    })
  }

### direct repeats ###
  observeEvent(input$strain, {
    create_direct_repeats_data(load_dataset(),
      input$strain,
      input$selected_segment,
      input$direct_repeats_flattened
    )
    output$direct_repeats_plot <- renderPlotly({
      create_direct_repeats_plot()
    })
  })
  observeEvent(input$selected_segment, {
    create_direct_repeats_data(load_dataset(),
      input$strain,
      input$selected_segment,
      input$direct_repeats_flattened
    )
    output$direct_repeats_plot <- renderPlotly({
      create_direct_repeats_plot()
    })
  })
  observeEvent(input$direct_repeats_flattened, {
    create_direct_repeats_data(load_dataset(),
      input$strain,
      input$selected_segment,
      input$direct_repeats_flattened
    )
    output$direct_repeats_plot <- renderPlotly({
      create_direct_repeats_plot()
    })
  })

### regression ###
  output$regression_plot <- renderPlot({
    create_regression_plot(
      load_dataset(),
      input$strain,
      input$regression_segments
    )
  })

}

###########
### APP ###
###########
# Main function that runs some prechecks and then builds the application
run_prechecks()
shinyApp(ui=ui, server=server)

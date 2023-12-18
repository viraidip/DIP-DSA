library(bslib)
library(comprehenr)
library(DT)
library(ggplot2)
library(hash)
library(jsonlite)
library(plotly)
library(plyr)
library(reticulate)
library(shiny)
library(shinydashboard)
library(shinyvalidate)
library(stringr)
library(tools)
library(ggvenn)
library(ComplexHeatmap)

library("Biostrings")

source("utils.R")

##########
### UI ###
##########
# Loads the sources for the UI of each tab.
# Each tab is saved in an individual file.
source('ui/new_dataset.R', local=TRUE)

source('ui/single_dataset.R', local=TRUE)
source('ui/multiple_datasets.R', local=TRUE)
source('ui/dataset_intersection.R', local=TRUE)

source('ui/about.R', local=TRUE)

ui <- bootstrapPage(
  dashboardPage(
    title="DIP-Deletion Site Analyzer",
    skin="green",
    dashboardHeader(title="DIP-DSA"),
    dashboardSidebar(
      sidebarMenu(
        id="sidebarmenu",
        menuItem("Add new dataset",
          tabName="new_dataset",
          icon=icon("database")
        ),
        hr(),
        menuItem("Single dataset",
          tabName="single_dataset",
          icon=icon("ruler-horizontal")
        ),
        menuItem("Multiple datasets",
          tabName="multiple_datasets",
          icon=icon("magnifying-glass-chart")
        ),
        menuItem("Dataset intersection",
          tabName="dataset_intersection",
          icon=icon("magnifying-glass")
        ),
        hr(),
        menuItem("About",
          tabName="about",
          icon=icon("circle-info")
        )
      )
    ),
    dashboardBody(
      tabItems(
        new_dataset_tab,
        
        single_dataset_tab,
        multiple_datasets_tab,
        dataset_intersection_tab,
        
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
source("server/new_dataset.R", local=TRUE)

source("server/single_dataset.R", local=TRUE)
source("server/multiple_datasets.R", local=TRUE)
source("server/dataset_intersection.R", local=TRUE)

source("server/about.R", local=TRUE)

server <- function(input, output, session) {
### new dataset ###
  observeEvent(input$link_to_about_tab, {
    updateTabItems(session, "sidebarmenu", "about")
  })

  observeEvent(input$dataset_submit, {
    # check if all fields are filled
    req(input$upload_strain, input$upload_dataset, input$upload_dataset_file)
    upload_strain <- format_strain_name(input$upload_strain)
    strain_path <- file.path(DATASETSPATH, upload_strain)

    # check if .csv file already exists and rename if so
    dataset_name <- input$upload_dataset
    file_path <- file.path(strain_path, paste(dataset_name, ".csv", sep=""))
    idx <- 0
    while (file.exists(file_path)) {
      idx <- idx + 1
      f_name <- paste(dataset_name, "_", idx, ".csv", sep="")
      file_path <- file.path(strain_path, f_name)
    }
    to_list <- list(file_path)
    from_list <- list(input$upload_dataset_file$datapath)
    move_files(from_list, to_list)

    create_random_data(upload_strain, dataset_name)
  })

  observeEvent(input$strain_submit, {
    # check if all fields are filled
    req(
      input$new_strain, input$upload_PB2_file, input$upload_PB1_file,
      input$upload_PA_file, input$upload_HA_file, input$upload_NP_file,
      input$upload_NA_file, input$upload_M_file, input$upload_NS_file
    )

    # move submitted files to right folder
    from_list <- list(
      input$upload_PB2_file$datapath, input$upload_PB1_file$datapath,
      input$upload_PA_file$datapath, input$upload_HA_file$datapath,
      input$upload_NP_file$datapath, input$upload_NA_file$datapath,
      input$upload_M_file$datapath, input$upload_NS_file$datapath
    )

    # check if strain exists create a folder if not
    upload_strain <- format_strain_name(input$new_strain)
    strain_path <- file.path(DATASETSPATH, upload_strain)
    if (dir.exists(strain_path)) {
      return()     
    } else {
      dir.create(strain_path)
    }
    fasta_path <- file.path(strain_path, "fastas")
    dir.create(fasta_path)
    
    # create list with paths on where to save the files and then move them
    to_list <- list()
    for (s in SEGMENTS) {
      f_name <- paste(s, ".fasta", sep="")
      to_list <- append(to_list, file.path(fasta_path, f_name))
    }
    move_files(from_list, to_list)
  })



### single dataset ###
  # Loading and dataset selection
  observe({
    path <- file.path(DATASETSPATH, format_strain_name(input$single_strain))
    dataset_names <- tools::file_path_sans_ext(list.files(path, pattern="csv"))
    updateSelectInput(session, "dataset", choices=dataset_names)
  })

  load_dataset <- reactive({
    path <- file.path(
      DATASETSPATH,
      format_strain_name(input$single_strain),
      paste(input$single_dataset, ".csv", sep="")
    )
    names <- c("Segment", "Start", "End", "NGS_read_count")
    cl <- c("character", "integer", "integer", "integer")
    if (file.exists(path)) {
      df <- read.csv(path, na.strings=c("NaN"), col.names=names, colClasses=cl)
    } else {
      df <- data.frame(
        "Segment"=character(),
        "Start"=integer(),
        "End"=integer(),
        "NGS_read_count"=integer()
      )
    }
    return(df)
  })

  observeEvent(input$link_to_single_data_point_tab, {
    updateTabItems(session, "sidebarmenu", "single_data_point")
  })

  # NGS counts
  output$ngs_distribution_plot <- renderPlotly(
    plot_ngs_distribution(
      load_dataset(),
      input$single_dataset
    )
  )

  # segment distribution
  output$segment_distribution_plot <- renderPlotly({
    create_segment_distribution_plot(
      load_dataset(),
      input$single_flattened,
      input$RCS_single
    )
  })

  # deletion shifts
  output$deletion_shift_plot <- renderPlotly({
    create_deletion_shift_plot(
      load_dataset(),
      input$single_flattened,
      input$RCS_single
    )
  })

  # lengths  
  output$lengths_plot <- renderPlotly({
    create_lengths_plot(
      load_dataset(),
      format_strain_name(input$single_strain),
      input$single_selected_segment,
      input$single_flattened,
      input$single_lengths_bins,
      input$RCS_single
    )
  })

  output$locations_plot <- renderPlotly({
    create_locations_plot(
      load_dataset(),
      format_strain_name(input$single_strain),
      input$single_selected_segment,
      input$single_flattened,
      input$RCS_single
    )
  })

  output$start_end_connection_plot <- renderPlotly({
    create_start_end_connection_plot(
      load_dataset(),
      input$single_dataset,
      format_strain_name(input$single_strain),
      input$single_selected_segment,
      input$RCS_single
    )
  })
  
  output$end_3_5_plot <- renderPlotly({
    create_end_3_5_plot(
      load_dataset(),
      format_strain_name(input$single_strain),
      input$single_selected_segment,
      input$RCS_single
    )
  })


  # Nucleotide enrichment
  observeEvent(
    eventExpr = {
      input$single_strain
      input$single_dataset
      input$single_selected_segment
      input$single_flattened
    },
    handlerExpr = {
      update_nuc_dist_plots()
    }
  )

  # function is called, when one of the inputs is changed (lines above)
  update_nuc_dist_plots <- function() {
    output$nuc_dist_start_A <- renderPlotly({
      create_nuc_dist_plot("Start", "A", input$single_selected_segment, input$single_dataset, input$RCS_single, format_strain_name(input$single_strain), input$single_flattened)
    })
    output$nuc_dist_start_C <- renderPlotly({
      create_nuc_dist_plot("Start", "C", input$single_selected_segment, input$single_dataset, input$RCS_single, format_strain_name(input$single_strain), input$single_flattened)
    })
    output$nuc_dist_start_G <- renderPlotly({
      create_nuc_dist_plot("Start", "G", input$single_selected_segment, input$single_dataset, input$RCS_single, format_strain_name(input$single_strain), input$single_flattened)
    })
    output$nuc_dist_start_U <- renderPlotly({
      create_nuc_dist_plot("Start", "U", input$single_selected_segment, input$single_dataset, input$RCS_single, format_strain_name(input$single_strain), input$single_flattened)
    })
  
    output$nuc_dist_end_A <- renderPlotly({
      create_nuc_dist_plot("End", "A", input$single_selected_segment, input$single_dataset, input$RCS_single, format_strain_name(input$single_strain), input$single_flattened)
    })
    output$nuc_dist_end_C <- renderPlotly({
      create_nuc_dist_plot("End", "C", input$single_selected_segment, input$single_dataset, input$RCS_single, format_strain_name(input$single_strain), input$single_flattened)
    })
    output$nuc_dist_end_G <- renderPlotly({
      create_nuc_dist_plot("End", "G", input$single_selected_segment, input$single_dataset, input$RCS_single, format_strain_name(input$single_strain), input$single_flattened)
    })
    output$nuc_dist_end_U <- renderPlotly({
      create_nuc_dist_plot("End", "U", input$single_selected_segment, input$single_dataset, input$RCS_single, format_strain_name(input$single_strain), input$single_flattened)
    })
  }

  # direct repeats
  observeEvent(
    eventExpr = {
      input$single_strain
      input$single_dataset
      input$single_selected_segment
      input$single_flattened
    },
    handlerExpr = {
      output$direct_repeats_plot <- renderPlotly({
        create_direct_repeats_plot(
          format_strain_name(input$single_strain),
          input$single_dataset,
          input$single_selected_segment,
          input$RCS_single,
          input$single_flattened
        )
      })
    }
  )



### multiple datasets





### dataset intersection
  output$dataset_venn <- renderPlot({
    df2 <- data.frame() # remove later
    plot_venn(
      load_dataset(),
      df2,
      input$strain,
      input$dataset,
      input$strain2,
      input$dataset2
    )
  })

  output$candidate_intersection_table <- renderDataTable({
    df2 <- data.frame() # remove later
    datatable(intersecting_candidates(load_dataset(), df2))
  }
  )

  output$intersecting_candidates_table <- renderDataTable(
    datatable(all_intersecting_candidates(
      input$selected_datasets,
      input$min_occurrences
      ),
      selection="single"
    )
  )
  
  output$overlap_matrix_plot <- renderPlotly({
    plot_overlap_matrix(
      input$selected_datasets
    )
  })

  output$upset_plot <- renderPlot({
    plot_upset_plot(
      input$selected_datasets
    )
  })










### about ###
  output$dataset_info_table <- renderTable({
    create_dataset_info_table()
  })

}


###########
### APP ###
###########
# Main function that runs some prechecks and then builds the application
run_prechecks()
shinyApp(ui=ui, server=server)

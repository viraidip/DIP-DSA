library(DT)
library(ggplot2)
library(hash)
library(jsonlite)
library(plotly)
library(plyr)
library(shiny)
library(shinydashboard)
library(shinyvalidate)
library(stringr)
library(dplyr)

# these two are from Bioconductor
library("Biostrings")
library("ComplexHeatmap")

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
    skin="red",
    dashboardHeader(title="DIP-DSA"),
    dashboardSidebar(
      sidebarMenu(
        id="sidebarmenu",
        style="position:fixed;width:220px;",
        menuItem("Add new dataset",
          tabName="new_dataset",
          icon=icon("file-circle-plus")
        ),
        hr(),
        menuItem("Single dataset",
          tabName="single_dataset",
          icon=icon("folder")
        ),
        menuItem("Multiple datasets",
          tabName="multiple_datasets",
          icon=icon("folder-tree")
        ),
        menuItem("Dataset intersection",
          tabName="dataset_intersection",
          icon=icon("magnifying-glass-chart")
        ),
        hr(),
        menuItem("About",
          tabName="about",
          icon=icon("circle-info")
        )
      )
    ),
    dashboardBody(
      tags$style(HTML("
        body {
          font-size: 16px;
        }"
      )),
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
    dataset_progress <- shiny::Progress$new()
    dataset_progress$set(message="Start upload", value=0)

    # check if all fields are filled
    req(input$upload_strain, input$upload_dataset, input$upload_dataset_file)
    upload_strain <- format_strain_name(input$upload_strain)
    strain_path <- file.path(DATASETSPATH, upload_strain)

    # check if .csv file already exists and rename if so
    dataset_progress$set(message="Check filename", value=0.1)
    dataset_name <- input$upload_dataset
    f_name <- dataset_name
    file_path <- file.path(strain_path, paste(f_name, ".csv", sep=""))
    idx <- 0
    while (file.exists(file_path)) {
      idx <- idx + 1
      f_name <- paste(dataset_name, "_", idx, sep="")
      file_path <- file.path(strain_path, paste(f_name, ".csv", sep=""))
    }
    from_list <- list(input$upload_dataset_file$datapath)
    to_list <- list(file_path)
    move_files(from_list, to_list)

    dataset_progress$set(message="Create random data", value=0.2)
    df <- read.csv(file_path)
    if (nrow(df) > 0) {
      create_random_data(upload_strain, f_name, dataset_progress)
    } else {
      to_list <- list(file.path(strain_path, paste(f_name, ".tsv", sep="")))
      move_files(from_list, to_list)
    }
    

    # add new options to select input
    updateSelectInput(
      session,
      inputId="single_strain",
      choices=gsub(
        "_",
        "/",
        list.dirs(DATASETSPATH, full.names=FALSE, recursive=FALSE)
      )
    )
    d_sets <- tools::file_path_sans_ext(list.files(strain_path, pattern="csv"))
    updateSelectInput(
      session,
      inputId="single_dataset",
      choices=d_sets
    )

    c <- list.files(DATASETSPATH, "csv$", full.names=FALSE, recursive=TRUE)
    updateSelectInput(
      session,
      inputId="multiple_datasets",
      choices=c,
      selected=c[1:2]
    )
    updateSelectInput(
      session,
      inputId="selected_datasets",
      choices=c,
      selected=c[1:2]
    )
    dataset_progress$close()
  })

  observeEvent(input$strain_submit, {
    strain_progress <- shiny::Progress$new()
    strain_progress$set(message="Start upload", value=0)
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

    strain_progress$set(message="Check strain name", value=0.4)
    # check if strain exists create a folder if not
    upload_strain <- format_strain_name(input$new_strain)
    strain_path <- file.path(DATASETSPATH, upload_strain)
    if (dir.exists(strain_path)) {
      strain_progress$close()
      return()     
    } else {
      dir.create(strain_path)
    }
    fasta_path <- file.path(strain_path, "fastas")
    dir.create(fasta_path)
    
    strain_progress$set(message="Upload FASTA files", value=0.4)
    # create list with paths on where to save the files and then move them
    to_list <- list()
    for (s in SEGMENTS) {
      f_name <- paste(s, ".fasta", sep="")
      to_list <- append(to_list, file.path(fasta_path, f_name))
    }
    move_files(from_list, to_list)

    c <- gsub("_","/",list.dirs(DATASETSPATH,full.names=FALSE,recursive=FALSE))
    updateSelectInput(
      session,
      inputId="upload_strain",
      choices=c
    )
    strain_progress$close()
  })


### single dataset ###
######################
  # update dataset selection when changing strain
  observe({
    path <- file.path(DATASETSPATH, format_strain_name(input$single_strain))
    dataset_names <- tools::file_path_sans_ext(list.files(path, pattern="csv"))
    updateSelectInput(session, "single_dataset", choices=dataset_names)
  })

  observeEvent(input$single_submit, {
    # initialize progess bar
    single_progress <- shiny::Progress$new()
    single_progress$set(message="Generating plots", value=0)

    # NGS counts
    output$ngs_distribution_plot <- renderPlotly({
      plot_ngs_distribution(
        format_strain_name(isolate(input$single_strain)),
        isolate(input$single_dataset),
        isolate(input$single_RSC),
        single_progress
      )
    })

    # frame shifts
    output$frame_shift_plot <- renderPlotly({
      plot_frame_shift(
        format_strain_name(isolate(input$single_strain)),
        isolate(input$single_dataset),
        isolate(input$single_flattened),
        isolate(input$single_RSC),
        single_progress
      )
    })

    # segment distribution
    output$segment_distribution_plot <- renderPlotly({
      plot_segment_distribution(
        format_strain_name(isolate(input$single_strain)),
        isolate(input$single_dataset),
        isolate(input$single_flattened),
        isolate(input$single_RSC),
        single_progress
      )
    })

    # lengths  
    output$lengths_plot <- renderPlotly({
      plot_lengths(
        format_strain_name(isolate(input$single_strain)),
        isolate(input$single_dataset),
        isolate(input$single_selected_segment),
        isolate(input$single_flattened),
        input$single_lengths_bins,
        isolate(input$single_RSC),
        single_progress
      )
    })
    
    # location
    output$locations_plot <- renderPlotly({
      plot_locations(
        format_strain_name(isolate(input$single_strain)),
        isolate(input$single_dataset),
        isolate(input$single_selected_segment),
        isolate(input$single_flattened),
        isolate(input$single_RSC),
        single_progress
      )
    })
    
    # compare 5' and 3'
    output$end_3_5_plot <- renderPlotly({
      plot_end_3_5(
        format_strain_name(isolate(input$single_strain)),
        isolate(input$single_dataset),
        isolate(input$single_selected_segment),
        isolate(input$single_RSC),
        single_progress
      )
    })

    # mapping start-end
    output$start_end_mapping_plot <- renderPlot({
      plot_start_end_mapping(
        format_strain_name(isolate(input$single_strain)),
        isolate(input$single_dataset),
        isolate(input$single_selected_segment),
        isolate(input$single_RSC),
        single_progress
      )
    })

    # direct repeats
    output$direct_repeats_plot <- renderPlotly({
      plot_direct_repeats(
        format_strain_name(isolate(input$single_strain)),
        isolate(input$single_dataset),
        isolate(input$single_selected_segment),
        isolate(input$single_RSC),
        isolate(input$single_flattened),
        single_progress
      )
    })

    # nucleotide enrichment
    output$nucleotide_enrichment_start_plot <- renderPlotly({
      plot_nucleotide_enrichment(
        format_strain_name(isolate(input$single_strain)),
        isolate(input$single_dataset),
        "Start",
        input$enrichment_nucleotide_start,
        isolate(input$single_selected_segment),
        isolate(input$single_RSC),
        isolate(input$single_flattened),
        single_progress
      )
    })
    output$nucleotide_enrichment_end_plot <- renderPlotly({
      plot_nucleotide_enrichment(
        format_strain_name(isolate(input$single_strain)),
        isolate(input$single_dataset),
        "End",
        input$enrichment_nucleotide_end,
        isolate(input$single_selected_segment),
        isolate(input$single_RSC),
        isolate(input$single_flattened),
        single_progress
      )
    })
  })


### multiple datasets ###
#########################
  observeEvent(input$multiple_submit, {
    # initialize progess bar
    multiple_progress <- shiny::Progress$new()
    multiple_progress$set(message="Generating plots", value=0)
    
    # NGS counts
    output$multiple_ngs_distribution_plot <- renderPlotly({
      plot_multiple_ngs_distribution(
        isolate(input$multiple_datasets),
        isolate(input$multiple_RSC),
        multiple_progress
      )
    })

    # segment distribution
    output$multiple_segment_distribution_plot <- renderPlotly({
      plot_multiple_segment_distribution(
        isolate(input$multiple_datasets),
        isolate(input$multiple_flattened),
        isolate(input$multiple_RSC),
        multiple_progress
      )
    })

    # deletion shifts
    output$multiple_deletion_shift_plot <- renderPlotly({
      plot_multiple_deletion_shift(
        isolate(input$multiple_datasets),
        isolate(input$multiple_flattened),
        isolate(input$multiple_RSC),
        multiple_progress
      )
    })

    # DVG length
    output$multiple_deletion_length_plot <- renderPlotly({
      plot_multiple_deletion_length(
        isolate(input$multiple_datasets),
        isolate(input$multiple_selected_segment),
        isolate(input$multiple_flattened),
        input$multiple_lengths_bins,
        isolate(input$multiple_RSC),
        multiple_progress
      )
    })

    # nucleotide enrichment
    output$multiple_nucleotide_enrichment_start_plot <- renderPlotly({
      plot_multiple_nucleotide_enrichment(
        isolate(input$multiple_datasets),
        isolate(input$multiple_selected_segment),
        "Start",
        isolate(input$multiple_flattened),
        input$multiple_enrichment_nucleotide_start,
        isolate(input$multiple_RSC),
        multiple_progress
      )
    })
    output$multiple_nucleotide_enrichment_end_plot <- renderPlotly({
      plot_multiple_nucleotide_enrichment(
        isolate(input$multiple_datasets),
        isolate(input$multiple_selected_segment),
        "End",
        isolate(input$multiple_flattened),
        input$multiple_enrichment_nucleotide_end,
        isolate(input$multiple_RSC),
        multiple_progress
      )
    })

    # direct repeats
    output$multiple_direct_repeats_plot <- renderPlotly({
      plot_multiple_direct_repeat(
        isolate(input$multiple_datasets),
        isolate(input$multiple_selected_segment),
        isolate(input$multiple_flattened),
        isolate(input$multiple_RSC),
        multiple_progress
      )
    })
  })


### dataset intersection ###
############################
  observeEvent(input$intersection_submit, {
    # initialize progess bar
    intersection_progress <- shiny::Progress$new()
    intersection_progress$set(message="Generating plots", value=0)

    # matrix with intersecting candidates
    output$overlap_matrix_plot <- renderPlotly({
      plot_overlap_matrix(
        isolate(input$selected_datasets),
        intersection_progress
      )
    })
    
    # barplot with candidates that occur in at least half of the datasets
    output$barplot_candidates_plot <- renderPlotly({
      plot_barplot_candidates(
        isolate(input$selected_datasets),
        intersection_progress
      )
    })

    # highest n ranked candidates
    output$highest_n_ranked_plot <- renderPlotly({
      plot_highest_n_ranked(
        isolate(input$selected_datasets),
        input$intersection_selected_segment,
        input$intersection_thresh,
        intersection_progress
      )
    })

    # table of highest n ranked candidates
    output$highest_n_ranked_table <- renderDataTable(
      datatable(
        get_highest_n_ranked_table(
          isolate(input$selected_datasets),
          input$intersection_selected_segment,
          input$intersection_thresh,
          intersection_progress
        )
      )
    )
  })



### about ###
  # No functions needed right now
}

###########
### APP ###
###########
# Main function that runs some prechecks and then builds the application
run_prechecks()
shinyApp(ui=ui, server=server)

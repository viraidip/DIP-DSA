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

library("Biostrings")

source("utils.R")
use_condaenv(CONDAENV)

##########
### UI ###
##########
# Loads the sources for the UI of each tab.
# Each tab is saved in an individual file.
source('ui/load_dataset.R', local=TRUE)
source('ui/dataset.R', local=TRUE)
source('ui/single_data_point.R', local=TRUE)
source('ui/lengths_locations.R', local=TRUE)
source('ui/nucleotide_distribution.R', local=TRUE)
source('ui/direct_repeats.R', local=TRUE)
source('ui/motif_search.R', local=TRUE)
source('ui/regression.R', local=TRUE)
source('ui/np_density.R', local=TRUE)
source('ui/candidate_prediction.R', local=TRUE)
source('ui/about.R', local=TRUE)

ui <- bootstrapPage(
  dashboardPage(
    title="DIP-Deletion Site Analyzer",
    skin="green",
    dashboardHeader(title="DIP-DSA"),
    dashboardSidebar(
      sidebarMenu(
        id="sidebarmenu",
        menuItem("(Up-)load dataset",
          tabName="load_dataset",
          icon=icon("database")
        ),
        menuItem("Dataset overview",
          tabName="dataset",
          icon=icon("table")
        ),
        menuItem("Inspect single data point",
          tabName="single_data_point",
          icon=icon("check")
        ),
        hr(),
        selectInput(
          inputId="selected_segment",
          label="Select segment",
          choices=SEGMENTS
        ),
        menuItem("Lengths and locations",
          tabName="lengths_locations",
          icon=icon("ruler-horizontal")
        ),
        menuItem("Nucleotide distribution",
          tabName="nucleotide_distribution",
          icon=icon("magnifying-glass-chart")
        ),
        menuItem("Direct repeats",
          tabName="direct_repeats",
          icon=icon("repeat")
        ),
        menuItem("Motif search",
          tabName="motif_search",
          icon=icon("magnifying-glass")
        ),
        menuItem("Nucleoprotein density",
          tabName="np_density",
          icon=icon("cubes-stacked")
        ),
        hr(),
        menuItem("Linear regression",
          tabName="regression",
          icon=icon("chart-line")
        ),
        menuItem("Candidate prediction",
          tabName="candidate_prediction",
          icon=icon("clipboard-check")
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
        load_dataset_tab,
        dataset_tab,
        single_data_point_tab,
        lengths_locations_tab,
        nucleotide_distribution_tab,
        direct_repeats_tab,
        motif_search_tab,
        regression_tab,
        np_density_tab,
        candidate_prediction_tab,
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
source("server/load_dataset.R", local=TRUE)
source("server/dataset.R", local=TRUE)
source("server/single_data_point.R", local=TRUE)
source("server/lengths_locations.R", local=TRUE)
source("server/nucleotide_distribution.R", local=TRUE)
source("server/direct_repeats.R", local=TRUE)
source("server/motif_search.R", local=TRUE)
source("server/regression.R", local=TRUE)
source("server/np_density.R", local=TRUE)
source("server/candidate_prediction.R", local=TRUE)
source("server/about.R", local=TRUE)

server <- function(input, output, session) {
### load/select dataset ###
  observe({
    path <- file.path(DATASETSPATH, format_strain_name(input$strain))
    dataset_names <- tools::file_path_sans_ext(list.files(path, pattern="csv"))
    updateSelectInput(session, "dataset", choices=dataset_names)
  })
  observe({
    path <- file.path(DATASETSPATH, format_strain_name(input$strain2))
    dataset_names <- tools::file_path_sans_ext(list.files(path, pattern="csv"))
    updateSelectInput(session, "dataset2", choices=dataset_names)
  })

  load_dataset <- reactive({
    path <- file.path(
      DATASETSPATH,
      format_strain_name(input$strain),
      paste(input$dataset, ".csv", sep="")
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

  observeEvent(input$link_to_about_tab, {
    updateTabItems(session, "sidebarmenu", "about")
  })

  observeEvent(input$dataset_submit, {
    # check if all fields are filled
    req(
      input$upload_strain, input$upload_dataset,
      input$upload_dataset_file,
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

    # check if strain exists create a folder if not
    upload_strain <- format_strain_name(input$upload_strain)
    strain_path <- file.path(DATASETSPATH, upload_strain)
    if (!dir.exists(strain_path)) {
      dir.create(strain_path)
      update_dataset <- FALSE
    } else {
      update_dataset <- TRUE
    }

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

    # check if fastas already exists and create folder if not
    fasta_path <- file.path(strain_path, "fastas")
    if (!dir.exists(fasta_path)) {
      dir.create(fasta_path)
    }

    # create list with paths on where to save the files and then move them
    for (s in SEGMENTS) {
      f_name <- paste(s, ".fasta", sep="")
      to_list <- append(to_list, file.path(fasta_path, f_name))
    }
    move_files(from_list, to_list)

    c <- gsub("_","/",list.dirs(DATASETSPATH,full.names=FALSE,recursive=FALSE))
    # select the new submitted dataset
    updateSelectInput(
      session,
      inputId="strain",
      choices=c,
      selected=input$upload_strain
    )
    if (update_dataset) {
      c <- tools::file_path_sans_ext(list.files(strain_path, pattern="csv"))
      updateSelectInput(
        session,
        inputId="dataset",
        choices=c,
        selected=dataset_name
      )
    }
  })


### dataset overview ###
  observeEvent(input$link_to_single_data_point_tab, {
    updateTabItems(session, "sidebarmenu", "single_data_point")
  })

  output$dataset_stats_info <- renderText(
    generate_stats_info(load_dataset()
    )
  )

  output$dataset_table <- renderDataTable(
    datatable(load_dataset(), selection="single")
  )

  output$dataset_venn <- renderPlot({
    df2 <- check_second_dataset(
      input$two_datasets,
      format_strain_name(input$strain2),
      input$dataset2
    )
    plot_venn(
      load_dataset(),
      df2,
      input$strain,
      input$dataset,
      input$strain2,
      input$dataset2
    )
  })


### single data point ###
  observeEvent(input$link_to_dataset_tab, {
    updateTabItems(session, "sidebarmenu", "dataset")
  })

  output$single_data_point_info <- renderText({
    create_single_data_point_info(
      load_dataset(),
      input$dataset_table_rows_selected,
      format_strain_name(input$strain)
    )
  })

  output$single_data_point_packaging_signal_info <- renderText({
    create_single_data_point_packaging_signal_info(
      load_dataset(),
      input$dataset_table_rows_selected,
      format_strain_name(input$strain)
    )
  })

  output$single_data_point_start_window <- renderPlot({
    plot_deletion_site_window(
      load_dataset(),
      input$dataset_table_rows_selected,
      format_strain_name(input$strain),
      "Start"
    )
  })

  output$single_data_point_end_window <- renderPlot({
    plot_deletion_site_window(
      load_dataset(),
      input$dataset_table_rows_selected,
      format_strain_name(input$strain),
      "End"
    )
  })


### lenghts and locations ###
  output$lengths_plot <- renderPlotly({
    df2 <- check_second_dataset(
      input$two_datasets,
      format_strain_name(input$strain2),
      input$dataset2
    )
    create_lengths_plot(
      load_dataset(),
      format_strain_name(input$strain),
      df2,
      input$strain2,
      input$selected_segment,
      input$lengths_flattened,
      input$lengths_bins
    )
  })

  output$locations_plot <- renderPlotly({
    df2 <- check_second_dataset(
      input$two_datasets,
      format_strain_name(input$strain2),
      input$dataset2
    )
    create_locations_plot(
      load_dataset(),
      df2,
      format_strain_name(input$strain),
      input$selected_segment,
      input$locations_flattened
    )
  })

  output$end_3_5_plot <- renderPlotly({
    create_end_3_5_plot(
      load_dataset(),
      format_strain_name(input$strain),
      input$selected_segment
    )
  })


### nucleotide distribution ###
  observeEvent(
    eventExpr = {
      input$strain
      input$dataset
      input$strain2
      input$dataset2
      input$two_datasets
      input$selected_segment
      input$nuc_dist_flattened
    },
    handlerExpr = {
      df2 <- check_second_dataset(
        input$two_datasets,
        format_strain_name(input$strain2),
        input$dataset2
      )
      create_nuc_dist_data(
        load_dataset(),
        format_strain_name(input$strain),
        df2,
        format_strain_name(input$strain2),
        input$selected_segment,
        input$nuc_dist_flattened
      )
      update_nuc_dist_plots()
    }
  )

  # function is called, when one of the inputs is changed (lines above)
  update_nuc_dist_plots <- function() {
    output$nuc_dist_start_A <- renderPlotly({
      create_nuc_dist_plot("Start", "A", input$selected_segment)
    })
    output$nuc_dist_start_C <- renderPlotly({
      create_nuc_dist_plot("Start", "C", input$selected_segment)
    })
    output$nuc_dist_start_G <- renderPlotly({
      create_nuc_dist_plot("Start", "G", input$selected_segment)
    })
    output$nuc_dist_start_U <- renderPlotly({
      create_nuc_dist_plot("Start", "U", input$selected_segment)
    })
  
    output$nuc_dist_end_A <- renderPlotly({
      create_nuc_dist_plot("End", "A", input$selected_segment)
    })
    output$nuc_dist_end_C <- renderPlotly({
      create_nuc_dist_plot("End", "C", input$selected_segment)
    })
    output$nuc_dist_end_G <- renderPlotly({
      create_nuc_dist_plot("End", "G", input$selected_segment)
    })
    output$nuc_dist_end_U <- renderPlotly({
      create_nuc_dist_plot("End", "U", input$selected_segment)
    })
  }


### direct repeats ###
  observeEvent(
    eventExpr = {
      input$strain
      input$dataset
      input$strain2
      input$dataset2
      input$two_datasets
      input$selected_segment
      input$direct_repeats_flattened
    },
    handlerExpr = {
      df2 <- check_second_dataset(
        input$two_datasets,
        format_strain_name(input$strain2),
        input$dataset2
      )
      create_direct_repeats_data(
        load_dataset(),
        format_strain_name(input$strain),
        df2,
        format_strain_name(input$strain2),
        input$selected_segment,
        input$direct_repeats_flattened
      )
      output$direct_repeats_plot <- renderPlotly({
        create_direct_repeats_plot(
          input$direct_repeats_correction,
          input$selected_segment
        )
      })
    }
  )

  output$direct_repeats_plot <- renderPlotly({
    create_direct_repeats_plot(input$direct_repeats_correction)
  })


### motif search ###
  output$motif_on_sequence <- renderPlotly({
    create_motif_on_sequence_plot(
      load_dataset(),
      format_strain_name(input$strain),
      input$selected_segment,
      input$motif,
      input$mismatch
    )
  })

  observeEvent(
    eventExpr = {
      input$strain
      input$dataset
      input$selected_segment
      input$motif
      input$mismatch
    },
    handlerExpr = {
      output$motif_table <- renderDataTable(
        create_motif_table(
          load_dataset(),
          format_strain_name(input$strain),
          input$selected_segment,
          input$motif,
          input$mismatch
        )
      )
    }
  )


### NP density ###
  observeEvent(input$selected_segment, {
    updateTextInput(
      session,
      inputId="np_areas",
      value=NP_MAP[[input$selected_segment]]
    )
  })

  observeEvent(
    eventExpr = {
      input$strain
      input$dataset
      input$np_areas
    },
    handlerExpr = {
      output$np_plot <- renderPlotly({
        create_np_plot(
          load_dataset(),
          format_strain_name(input$strain),
          input$selected_segment,
          input$np_areas
        )
      })
    }
  )

  output$np_ratios_info <- renderText({
    create_np_ratios_info(
      load_dataset(),
      format_strain_name(input$strain),
      input$selected_segment,
      input$np_areas
    )
  })


### regression ###
  output$regression_plot <- renderPlot({
    create_regression_plot(
      load_dataset(),
      format_strain_name(input$strain),
      input$regression_segments
    )
  })


### candidate prediction ###
  observe({
    max <- get_seq_len(format_strain_name(input$prediction_strain),
      input$prediction_segment)
    iv <- InputValidator$new()
    iv$add_rule("prediction_start", sv_between(1, max))
    iv$add_rule("prediction_end", sv_between(1, max))
    iv$enable()
  })
  observeEvent(
    eventExpr = {
      input$run_prediction
    },
    handlerExpr = {
      label <- run_prediction(
        input$prediction_start,
        input$prediction_end,
        input$prediction_strain,
        input$prediction_segment,
        input$classifier
      )
      segment <- input$prediction_segment
      start <- input$prediction_start
      end <- input$prediction_end
      output$prediction_results <- renderText({
        show_prediction_results(
          label, segment, start, end
        )
      })
    }
  )


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

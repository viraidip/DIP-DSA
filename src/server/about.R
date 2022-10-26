
create_dataset_info_table <- function() {
  path <- file.path(DATAPATH, "datasets_metadata.csv")
  return(read.csv(path))
}


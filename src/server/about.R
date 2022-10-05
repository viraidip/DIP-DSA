
create_dataset_info_table <- function() {
  path <- file.path(DATAPATH, "datasets_metadata.csv")
  read.csv(path)
}

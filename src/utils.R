
SEGMENTS <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")

DATAPATH <- file.path("..", "data")
DATASETSPATH <- file.path(DATAPATH, "datasets")
FASTAPATH <- file.path(DATAPATH, "strain_segment_fastas")
TEMPPATH <- file.path(DATAPATH, "temp")

COLOR_MAP <- hash(A="blue", C="yellow", G="green", U="red")
NUC_MAP <- hash(A="Adenin", C="Cytosin", G="Guanin", U="Uracil")

run_prechecks <- function() {
  if (!file.exists(TEMPPATH)) {
    dir.create(TEMPPATH)
  }
}

get_seq <- function(strain, segment) {
  path <- file.path(FASTAPATH, strain, paste(segment, ".fasta", sep=""))
  if (!file.exists(path)) {
    return(NULL)
  }
  fasta <- readDNAStringSet(path)
  return(RNAString(fasta[[1]]))
}

get_seq_len <- function(strain, segment) {
  return(length(get_seq(strain, segment)))
}

get_stat_symbol <- function(p) {
  if (p < 0.0001) {
    return("****")
  }
  else if (p < 0.001) {
    return("***")
  }
  else if (p < 0.01) {
    return("**")
  }
  else if (p < 0.05) {
    return("*")
  }
  else {
    return("")
  }
}

packaging_signal_data_exists <- function(strain) {
  path <- file.path(DATAPATH, "packaging_signal", paste(strain, ".json", sep=""))
  return(file.exists(path))
}

load_packaging_signal_data <- function(strain) {
  path <- file.path(DATAPATH, "packaging_signal", paste(strain, ".json", sep=""))
  data <- read_json(path)
  return(data)
}


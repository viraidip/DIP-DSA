### DEFINING STATIC VALUES ###
SEGMENTS <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")

DATAPATH <- file.path("..", "data")
DATASETSPATH <- file.path(DATAPATH, "datasets")
FASTAPATH <- file.path(DATAPATH, "strain_segment_fastas")
TEMPPATH <- file.path(DATAPATH, "temp")

COLOR_MAP <- hash(A="blue", C="yellow", G="green", U="red")
NUC_MAP <- hash(A="Adenin", C="Cytosin", G="Guanin", U="Uracil")

### DEFINING FUNCTIONS ###
run_prechecks <- function() {
  if (!file.exists(TEMPPATH)) {
    dir.create(TEMPPATH)
  }
}

get_seq <- function(strain, segment) {
  path <- file.path(DATASETSPATH, strain, "fastas", paste(segment, ".fasta", sep=""))
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
  if (is.nan(p)) {
    return("")
  } else if (p < 0.0001) {
    return("****")
  } else if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("")
  }
}

packaging_signal_data_exists <- function(strain) {
  path <- file.path(DATASETSPATH, strain, "packaging_signal.json")
  return(file.exists(path))
}

load_packaging_signal_data <- function(strain) {
  path <- file.path(DATASETSPATH, strain, "packaging_signal.json")
  data <- read_json(path)
  return(data)
}

# reformat name of strain folder
# needs to be called everytime the strain folder is accessed
# folder is e.g. A_California_07_2009
# but it is displayed and can be selected as A/California/07/2009
format_strain_name <- function(strain) {
  return(gsub(pattern="/", replacement="_", x=strain))
}

check_second_dataset <- function(include, strain, dataset) {
  if (include == "Yes") {
    path <- file.path(DATASETSPATH, strain, paste(dataset, ".csv", sep=""))
    names <- c("Segment", "Start", "End", "NGS_read_count")
    classes <- c("character", "integer", "integer", "integer")
    if (file.exists(path)) {
      df <- read.csv(path, na.strings=c("NaN"), col.names=names, colClasses=classes)
    }
  } else {
    df <- data.frame()
  }
  return (df)
}


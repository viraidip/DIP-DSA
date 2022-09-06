library("Biostrings")

SEGMENTS <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")

DATASETSPATH <- file.path("..", "data", "datasets")
FASTAPATH <- file.path("..", "data", "strain_segment_fastas")
TEMPPATH <- file.path("..", "data", "temp")

run_prechecks <- function() {
  if (!file.exists(TEMPPATH)) {
    dir.create(TEMPPATH)
  }
}

get_seq <- function(strain, segment) {
  path <- file.path(FASTAPATH, strain, paste(segment, ".fasta", sep=""))
  fasta <- readDNAStringSet(path)
  return(RNAString(fasta[[1]]))
}

get_seq_len <- function(strain, segment) {
  return(length(get_seq(strain, segment)))
}


library("Biostrings")

SEGMENTS <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")

get_seq <- function(strain, segment) {
  path <- file.path("..", "data", "strain_segment_fastas", strain, paste(segment, ".fasta", sep=""))
  fasta <- readDNAStringSet(path)
  return(fasta[[1]])
}

get_seq_len <- function(strain, segment) {
  path <- file.path("..", "data", "strain_segment_fastas", strain, paste(segment, ".fasta", sep=""))
  fasta <- readDNAStringSet(path)
  return(width(fasta))
}

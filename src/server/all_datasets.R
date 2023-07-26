
load_dataset <- function(fname) {
  path <- file.path(DATASETSPATH, fname)
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
}

all_intersecting_candidates <- function(datasets, thresh) {
  candidates <- c()
  for (p in datasets) {
    df <- load_dataset(p)
    cand <- paste(df$Segment, df$Start, df$End, sep="_")
    candidates <- append(candidates, cand)
  }
  counts <- table(candidates)
  selected <- counts[counts >= thresh]
  return(data.frame(selected))
}

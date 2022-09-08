
create_single_datapoint_info <- function(df, row, strain) {
  if (is.null(row)) {
    paste("No data point selected.")
  }
  else {
    segment <- df[row, "Segment"]
    start <- df[row, "Start"]
    end <- df[row, "End"]
    full_length <- get_seq_len(strain, segment)
    return(
      paste(
        paste("Selected datapoint belongs to segment", segment),
        paste("Start position is", start),
        paste("End position is", end),
        paste("Deletion length is" , end-start),
        paste("DIP sequence is of length", start+(full_length-end+1))
        ,sep="\n"
      )
    )
  }
}


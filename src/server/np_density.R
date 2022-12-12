
create_motif_on_sequence_plot <- function(df, strain, segment, motif) {
  df <- format_dataframe_locations(df, segment, "flattened")
  if (nrow(df) == 0) {
    return()
  }

  # search for motif on sequence
  if (nchar(motif) > 0) {
    sequence <- get_seq(strain, segment)
    matches <- matchPattern(motif, sequence)
  } else {
    matches <- c()
  }

  # create plot
  p <- ggplot(df, aes(x=Position, y=NGS_read_count, fill=Class)) +
    geom_bar(stat="identity") +
    xlab("Nucleotide position on segment") +
    ylab("NGS read count")

  # only draw if there are less than 100 matches
  if (length(matches) > 0 && length(matches) < 100) {
    c <- "black"
    for (i in 1:length(matches)) {
      m <- as(matches[i], "IRanges")
      xmin <- start(m)
      xmax <- end(m)
      p <- p + geom_rect(xmin=xmin, xmax=xmax, ymin=0, ymax=1, col=c, fill=c)
    }
  }

  ggplotly(p)
}

create_motif_table <- function(df, strain, segment, motif) {
  df <- format_dataframe_locations(df, segment, "flattened")
  if (nrow(df) > 0 && nchar(motif) > 0) {
    sequence <- get_seq(strain, segment)
    matches <- matchPattern(motif, sequence)
    m <- data.frame(as(matches, "IRanges"))
  } else {
    m <- data.frame()
  }
  return(m)
}



find_matches <- function (strain, segment, motif) {
  motif <- reformat_motif(motif)
  if (nchar(motif) > 0) {
    sequence <- get_seq(strain, segment)
    matches <- matchPattern(motif, sequence)
  } else {
    matches <- c()
  }

  return(matches)
}

create_motif_on_sequence_plot <- function(df, strain, segment, motif) {
  df <- format_dataframe_locations(df, segment, "flattened")
  if (nrow(df) == 0) {
    return()
  }

  # search for motif on sequence
  matches <- find_matches(strain, segment, motif)

  # create plot
  p <- ggplot(df, aes(x=Position, y=NGS_read_count, fill=Class)) +
    geom_bar(stat="identity", position="dodge", width=1) +
    xlab("Nucleotide position on segment") +
    ylab("NGS read count")

  # only draw if there are less than 100 matches
  if (length(matches) > 0 && length(matches) < 100) {
    matches <- data.frame(matches)
    matches["f"] <- "found motif"
    p <- p + geom_rect(data=matches,
      aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=f),
      inherit.aes=FALSE
    )
  }

  ggplotly(p)
}

create_motif_table <- function(df, strain, segment, motif) {
  df <- format_dataframe_locations(df, segment, "flattened")

  m <- data.frame()
  if (nrow(df) > 0) {
    matches <- find_matches(strain, segment, motif)
    if (length(matches) != 0) {
      m <- data.frame(as(matches, "IRanges"))
    }
  }
  return(m)
}


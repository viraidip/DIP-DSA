
find_matches <- function (strain, segment, motif, mismatch) {
  motif <- reformat_motif(motif)
  if (nchar(motif) > 0) {
    sequence <- get_seq(strain, segment)
    matches <- matchPattern(motif, sequence, max.mismatch=mismatch)
  } else {
    matches <- c()
  }

  return(matches)
}

create_motif_on_sequence_plot <- function(df,strain,segment,motif,mismatch) {
  df <- format_dataframe_locations(df, segment, "flattened")
  if (nrow(df) == 0) {
    return()
  }

  # search for motif on sequence
  matches <- find_matches(strain, segment, motif, mismatch)

  # create plot
  p <- ggplot(df, aes(x=Position, y=NGS_read_count, fill=Class)) +
    geom_bar(stat="identity", position="dodge", width=1) +
    xlab("Nucleotide position on segment") +
    ylab("NGS read count") +
    xlim(0, get_seq_len(strain, segment)) +
    ggtitle(paste("Matches of motif '",
      reformat_motif(motif),
      "' together with deletion sites for segment ",
      segment,
      sep="")
    ) +
    theme(plot.title = element_text(size=20))

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

create_motif_table <- function(df, strain, segment, motif, mismatch) {
  df <- format_dataframe_locations(df, segment, "flattened")

  m <- data.frame()
  if (nrow(df) > 0) {
    matches <- find_matches(strain, segment, motif, mismatch)
    if (length(matches) != 0) {
      if (mismatch == 0) {
        m <- data.frame(as(matches, "IRanges"))
      } else { # include sequence of matches
        m <- data.frame(matches)
        colnames(m)[4] <- "sequence"
      }
      m <- subset(m, select=-c(width))
    }
  }
  return(m)
}


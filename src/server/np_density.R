
create_np_plot <- function(df, strain, segment, areas) {
  # reformat data set
  df <- format_dataframe_locations(df, segment, "flattened")

  # reformat np areas
  a_df <- reformat_np_areas(areas)

  p <- ggplot(df, aes(x=Position, y=NGS_read_count, fill=Class)) +
    geom_bar(stat="identity") +
    xlab("Nucleotide position on segment") +
    ylab("NGS read count")

  if (nrow(a_df) > 0) {
    c <- "black"
    for (i in 1:nrow(a_df)) {
      xmin <- a_df[i, ]$start
      xmax <- a_df[i, ]$end
      p <- p + geom_rect(xmin=xmin, xmax=xmax, ymin=0, ymax=1, fill=c, alpha=0.3)
    }
  }

  ggplotly(p)
}

create_np_bar_plot <- function(df, strain, areas) {
  # reformat np areas


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
      p <- p + geom_rect(xmin=xmin, xmax=xmax, ymin=0, ymax=1, fill=c)
    }
  }

  ggplotly(p)
}


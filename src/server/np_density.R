
create_np_plot <- function(df, strain, segment, areas) {
  # reformat data
  df <- format_dataframe_locations(df, segment, "flattened")
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

in_high_np_area <- function(row, a_df) {
  p <- as.integer(row["Position"])
  for (i in 1:nrow(a_df)) {
    s <- a_df[i, "start"]
    e <- a_df[i, "end"]
    if ((p >= s) && (p <= e)) {
      return("high")
    }
  }
  return("low")
}

create_np_bar_plot <- function(df, strain, segment, areas) {
  # reformat data
  obs_df <- format_dataframe_locations(df, segment, "flattened")
  a_df <- reformat_np_areas(areas)

  seq <- get_seq(strain, segment)
  n_samples <- nrow(df) * 5
  sam_df <- create_direct_repeat_sampling_data(df, n_samples, seq)
  sam_df["Segment"] <- rep(segment, n_samples)
  sam_df <- format_dataframe_locations(sam_df, segment, "flattened")

  obs_df["np_area"] <- apply(obs_df, 1, in_high_np_area, a_df=a_df)
  sam_df["np_area"] <- apply(sam_df, 1, in_high_np_area, a_df=a_df)

  # create plot
  p <- ggplot(df, aes(x=Position, y=NGS_read_count, fill=Class)) +
    geom_bar(stat="identity") +
    xlab("Nucleotide position on segment") +
    ylab("NGS read count")


  ggplotly(p)
}



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
    if (p >= s && p <= e) {
      return("high")
    }
  }
  return("low")
}

count_area_occurrence <- function(df, a_df, label) {
  df["NP_area"] <- apply(df, 1, in_high_np_area, a_df=a_df)
  table_df <- data.frame(table(df$NP_area))
  # dividing high/low
  final_df <- data.frame(Var1="Ratio", Freq=table_df[1, 2]/table_df[2, 2])
  final_df["Class"] <- label
  return(final_df)
}

create_np_bar_plot <- function(df, strain, segment, areas) {
  # reformat data
  obs_df <- format_dataframe_locations(df, segment, "flattened")

  # create sampling data
  seq <- get_seq(strain, segment)
  n_samples <- nrow(df) * 5
  sam_df <- create_direct_repeat_sampling_data(df, n_samples, seq)
  sam_df["Segment"] <- rep(segment, n_samples)
  sam_df <- format_dataframe_locations(sam_df, segment, "flattened")

  # reformat np areas, if they are valid areas provided create a plot
  a_df <- reformat_np_areas(areas)
  if (nrow(a_df) == 0) {
    p <- ggplot()
  } else {
    # count occurrences inside and outside of NP areas
    o_df <- count_area_occurrence(obs_df, a_df, "observed")
    s_df <- count_area_occurrence(sam_df, a_df, "expected")
    plot_df <- rbind(o_df, s_df)

    # TODO: add statistical comparision

    # create plot
    p <- ggplot(data=plot_df, aes(x=Var1, y=Freq, fill=Class)) +
      geom_bar(stat="identity", position=position_dodge()) +
      xlab("Nucleotide position on segment") +
      ylab("NGS read count")
  }

  ggplotly(p)
}


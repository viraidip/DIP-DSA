create_locations_plot <- function(df, segment, flattened) {
  # slice df by segment, reformat and bind on position and NGS count
  df <- df[df$Segment == segment, ]
  start_df <- df[, c("Start", "NGS_read_count")]
  start_df["Class"] <- "Start"
  colnames(start_df)[colnames(start_df) == "Start"] <- "Position"
  end_df <- df[, c("End", "NGS_read_count")]
  end_df["Class"] <- "End"
  colnames(end_df)[colnames(end_df) == "End"] <- "Position"
  df <- rbind(start_df, end_df)

  # reformat if data should be displayed flattened
  if (flattened == "flattened") {
    df["NGS_read_count"] <- 1
  }

  p <- ggplot(df, aes(x=Position, y=NGS_read_count, fill=Class)) +
    geom_bar(stat="identity")
  ggplotly(p)
}

create_lengths_plot <- function(df, segment, strain, flattened, n_bins) {
  # slice df by segment, reformat and bind on position and NGS count
  df <- df[df$Segment == segment, ]

  seq_len <- get_seq_len(strain, segment)
  df["Length"] <- df["Start"] + (seq_len - df["End"] + 1)
  # multiply each column by NGS count if data is unflattened
  if (flattened != "flattened") {
    df <- data.frame(lapply(df, rep, df$NGS_read_count))
  }

  p <- ggplot(df, aes(x=Length)) +
    geom_histogram(binwidth=n_bins)
  ggplotly(p)
}

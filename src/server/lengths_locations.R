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

  ggplot(df, aes(x=Position, y=NGS_read_count, fill=Class)) +
    geom_bar(stat="identity")
}


create_lengths_plot <- function(df, segment, flattened, n_bins) {
  # slice df by segment, reformat and bind on position and NGS count
  df <- df[df$Segment == segment, ]
# TODO: here I need to include the full length segment length
  df["Length"] <- df["Start"] - df["End"]
  # reformat if data should be displayed flattened
  if (flattened == "flattened") {
    df["NGS_read_count"] <- 1
  }

  ggplot(df, aes(x=Length)) +
    geom_histogram(binwidth=n_bins)

}

create_lengths_plot <- function(df, segment, flattened) {
  # slice df by segment
  df <- df[df$segment == segment, ]
  # create barplot
  print(df)
  ggplot(df, aes(x=Start)) +
    geom_histogram(binwidth=1)
}

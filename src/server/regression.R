
create_regression_plot <- function(df, strain, segments) {
  df <- df[df$Segment %in% segments,]
  if (nrow(df) == 0) {
    return(0)
  }

  # reformat df, segment length x NGS_read_count
  data <- aggregate(list(relative_count=df$NGS_read_count), by=list(Segment=df$Segment), FUN=sum)
  data$relative_count <- data$relative_count / sum(data$relative_count)

  segment_lengths <- to_vec(for (seg in data$Segment) get_seq_len(strain, seg))
  data$segment_length <- segment_lengths
  regression_data <- data

  # do linear regression
  regression <- lm(relative_count ~ segment_length, data=regression_data)
  r_squared <- summary(regression)$r.squared
  slope <- regression$coefficients["segment_length"]
  intercept <- regression$coefficients["(Intercept)"]

  # get expected value by segment length
  expected_data <- data.frame(data)
  expected_data$relative_count <- segment_lengths/sum(segment_lengths)
  data$group <- rep("observed", nrow(data))
  expected_data$group <- rep("expected", nrow(expected_data))
  data <- rbind(data, expected_data)

  label <- paste("y =", format(intercept, digits=2), "+ x *", format(slope, digits=2), "\nR²:", format(r_squared, digits=2))

  # create the scatterplot
  ggplot(data, aes(x=segment_length, y=relative_count, color=group)) +
    geom_point() +
    geom_abline(intercept=intercept, slope=slope, show.legend=TRUE) +
    annotate(geom="text", x=500, y= 0.2, label=label) +
    xlim(0, max(data$segment_length)) +
    ylim(0, max(data$relative_count))

}

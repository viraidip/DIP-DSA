
create_regression_plot <- function(df, strain, segments) {
  df <- df[df$Segment %in% segments,]
  if (nrow(df) == 0) {
    return()
  }

  # reformat df, segment length x NGS_read_count
  data <- aggregate(list(relative_count=df$NGS_read_count),
    by=list(Segment=df$Segment),
    FUN=sum
  )
  data$relative_count <- data$relative_count / sum(data$relative_count)

  segment_lengths <- to_vec(for (seg in data$Segment) get_seq_len(strain, seg))
  data$segment_length <- segment_lengths
  regression_data <- data

  # do linear regression
  regression <- lm(relative_count ~ segment_length, data=regression_data)
  r_squared <- format(summary(regression)$r.squared, digits=2)
  slope <- regression$coefficients["segment_length"]
  intercept <- regression$coefficients["(Intercept)"]

  # get expected value by segment length
  expected_data <- data.frame(data)
  expected_data$relative_count <- segment_lengths/sum(segment_lengths)
  data$group <- rep("observed", nrow(data))
  expected_data$group <- rep("expected", nrow(expected_data))
  data <- rbind(data, expected_data)

  m <- format(slope, digits=2)
  c <- format(intercept, digits=2)
  func <- paste("f(x) =", m, "* x +", c)
  func_label <- paste(func, "\nR²:", r_squared)

  intersection <- -intercept/slope
  inter_label <- paste("(", format(intersection, digits=0), ", 0)", sep="")

  # create the scatterplot
  ggplot(
    data,
    aes(x=segment_length, y=relative_count, color=group, label=Segment)
  ) +
    geom_point() +
    geom_text(hjust=0, vjust=0, check_overlap=TRUE) +
    geom_abline(intercept=intercept, slope=slope, show.legend=TRUE) +
    annotate(geom="text", x=500, y=0.2, label=func_label) +
    geom_point(aes(x=intersection, y=0)) +
    annotate(geom="text", x=intersection, y=0, label=inter_label, hjust=0) +
    xlim(0, max(data$segment_length)) +
    ylim(0, max(data$relative_count)) +
    xlab("Length of segment") +
    ylab("Relative occurrence")
}


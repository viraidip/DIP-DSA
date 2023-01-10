
generate_stats_info <- function(df) {
  mean <- round(mean(df$NGS_read_count), 2)
  median <- median(df$NGS_read_count)
  max <- max(df$NGS_read_count)
  min <- min(df$NGS_read_count)
  return(
    paste(
      paste("Mean:\t", mean),
      paste("Median:\t", median),
      paste("Min:\t", min),
      paste("Max:\t", max),
      sep="\n"
    )
  )
}

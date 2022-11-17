
generate_stats_info <- function(df) {
  return(
    paste(
      paste("Mean:", mean(df$NGS_read_count)),
      paste("Median:", median(df$NGS_read_count)),
      paste("Max:", max(df$NGS_read_count)),
      paste("Min:", min(df$NGS_read_count)),
      sep="\n"
    )
  )
}

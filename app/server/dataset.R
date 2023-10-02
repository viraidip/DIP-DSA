
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

plot_ngs_distribution <- function(df1, d1) {
  p <- plot_ly(data = df1, y = ~NGS_read_count, type = "box", name=d1) %>%
    layout(
      title = "Boxplot of NGS read count (Log Scale)",
      yaxis = list(type = "log")
    )

  p
}

plot_venn <- function(df1, df2, s1, d1, s2, d2) {
  error_text <- "Select a second dataset to show plot."
  shiny::validate(need((nrow(df1 != 0) & (nrow(df2 != 0))), error_text))

  x <- list(
    A=paste(df1$Segment, df1$Start, df1$End, sep="_"),
    B=paste(df2$Segment, df2$Start, df2$End, sep="_")
  )
  names <- c(paste(s1, d1, sep=", "), paste(s2, d2, sep=", "))
  names(x) <- names

  p <- ggvenn(x, names, text_size=8) +
    ggtitle("Overlap of two selected datasets") + 
    theme(plot.title = element_text(size=20))

  p
}

intersecting_candidates <- function(df1, df2) {
  if (nrow(df1) == 0 || nrow(df2) == 0) {
    return(data.frame())
  }

  df1["DI"] <- paste(df1$Segment, df1$Start, df1$End, sep="_")
  df2["DI"] <- paste(df2$Segment, df2$Start, df2$End, sep="_")

  merged_df <- merge(df1, df2, by="DI", suffixes=c("_df1", "_df2"))
  final_dataset <- merged_df[, c("Segment_df1", "Start_df1", "End_df1", "NGS_read_count_df1", "NGS_read_count_df2")]
  colnames(final_dataset) <- c("Segment", "Start", "End", "NGS_read_count_dataset1", "NGS_read_count_dataset2")

  return(final_dataset)
}

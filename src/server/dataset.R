
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

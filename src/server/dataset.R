
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
  grid.newpage()
  if (nrow(df2) == 0) {
    grid.draw(x=NULL)
  } else {
  x <- list(
    A=paste(df1$Segment, df1$Start, df1$End, sep="_"),
    B=paste(df2$Segment, df2$Start, df2$End, sep="_")
  )
  venn_object <- venn.diagram(
    x,
    filename=NULL,
    disable.logging=TRUE,
    fill=c("#E69F00", "#56B4E9"),
    category.names=c(paste(s1, d1, sep="\n"), paste(s2, d2, sep="\n")),
    cex=1.5,
    cat.cex=1.5,
    cat.dist=c(-0.05,-0.05)
  )
  grid.draw(venn_object)
  }
}


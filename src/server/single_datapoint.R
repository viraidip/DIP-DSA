
create_single_datapoint_info <- function(df, row, strain) {
  if (is.null(row)) {
    paste("No data point selected.")
  }
  else {
    segment <- df[row, "Segment"]
    start <- df[row, "Start"]
    end <- df[row, "End"]
    ngs_read_count <- df[row, "NGS_read_count"]
    full_length <- get_seq_len(strain, segment)
    return(
      paste(
        paste("Selected datapoint belongs to segment", segment),
        paste("Start position is", start),
        paste("End position is", end),
        paste("Deletion length is" , end-start),
        paste("DIP sequence is of length", start+(full_length-end+1)),
        paste("Number of reads in sample (NGS read count) is", ngs_read_count)
        ,sep="\n"
      )
    )
  }
}

get_color <- function(nuc) {
  COLOR_MAP[[nuc]]
}

plot_deletion_site_window <- function(df, row, strain, pos) {
  segment <- df[row, "Segment"]
  # adjust by 1 if End site is given, to have indexing right afterwards
  position_int <- ifelse(pos == "End", df[row, pos]-1, df[row, pos])
  sequence <- get_seq(strain, segment)
  
  window <- subseq(sequence, start=position_int-4, end=position_int+5)
  nucleotide_list <- unlist(strsplit(as.character(window), split=""))
  color <- unlist(lapply(nucleotide_list, get_color))
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

  xleft <- ifelse(pos == "End", -1, 5.5)
  xright <- ifelse(pos == "End", 5.5, 12)

  if (pos == "End") {
    nucleotide_indexing <- c("...", paste("...", get_seq_len(strain, segment)))
  } else {
    nucleotide_indexing <- c("0 ...", "...")
  }

  plot(NULL, axes=FALSE, xlab="", ylab="", xlim=c(-1, 12), ylim=c(0, 2))
  # create boxes
  rect(xleft=-1, xright=12, ybottom=0, ytop=2)
  rect(xleft=xleft, xright=xright, ybottom=0, ytop=2, col="gray")
  # write single nucleotides
  text(x=x, y=1, label=nucleotide_list, col=color, cex=3)
  text(x=c(0, 11), y=1, label="...", cex=3)
  # generate descriptions
  text(x=x, y=0.2, label=x+position_int-5)
  text(x=c(0, 11), y=0.2, label=nucleotide_indexing)
  text(x=ifelse(pos == "End", 2, 9), y=1.9, "deleted sequence")
  text(x=ifelse(pos == "End", 9, 2), y=1.9, "remaining sequence")

}


library(hash)

counting_routine <- function(l, window, letter, ngs_read_count) {
  count_indices <- unlist(gregexpr(letter, window))
  if (count_indices[1] != -1){
    for (i in count_indices) {
      l[[i]] <- l[[i]] + ngs_read_count
    }
  }
  return (l)
}

count_nuc_dist <- function(seq, positions, ngs_read_counts) {
  A <- integer(10)
  C <- integer(10)
  G <- integer(10)
  U <- integer(10)
  for (i in 1:length(positions)) {
    p <- positions[[i]]
    ngs_read_count <- ngs_read_counts[[i]]
    window <- subseq(seq, start=p-4, end=p+5)
    A <- counting_routine(A, window, "A", ngs_read_count)
    C <- counting_routine(C, window, "C", ngs_read_count)
    G <- counting_routine(G, window, "G", ngs_read_count)
    U <- counting_routine(U, window, "U", ngs_read_count)
  }
  rel_occurrence <- c(A, C, G, U) / sum(ngs_read_counts)
  position <- c(rep(seq(1, 10), 4))
  nucleotide <- c(rep("A", 10), rep("C", 10), rep("G", 10), rep("U", 10))
  return(data.frame(rel_occurrence, position, nucleotide))
}

create_nuc_dist_plot <- function(df, strain, segment, pos, flattened, nuc) {
  # slice dataset
  df <- df[df$Segment == segment,]
  positions <- df[, pos]
  ngs_read_counts <- df[, "NGS_read_count"]
  if (flattened == "flattened") {
    ngs_read_counts[ngs_read_counts != 1] <- 1
  }

  # load sequence
  sequence <- get_seq(strain, segment)

  # count nuc dist around deletion site
  count_df <- count_nuc_dist(sequence, positions, ngs_read_counts)
  count_df <- count_df[count_df$nucleotide == nuc,]

  color <- hash()
  color[["A"]] <- "blue"
  color[["C"]] <- "green"
  color[["G"]] <- "yellow"
  color[["U"]] <- "red"

  # create a barplot
  ggplot(data=count_df, aes(x=position, y=rel_occurrence, fill=nucleotide)) +
    geom_bar(stat="identity", fill=color[[nuc]], position=position_dodge()) +
    scale_x_continuous(
      breaks=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
      labels=c("5", "4", "3", "2", "1", "-1", "-2", "-3", "-4", "-5")
    )

}

counting_routine <- function(l, window, letter) {
  count_indices <- unlist(gregexpr(letter, window))
  if (count_indices[1] != -1){
    for (i in count_indices) {
      l[[i]] <- l[[i]] + 1
    }
  }
  return (l)
}

count_nuc_dist <- function(seq, positions) {
  A <- integer(10)
  C <- integer(10)
  G <- integer(10)
  U <- integer(10)
  for (p in positions) {
    window <- subseq(seq, start=p-5, end=p+4)
    A <- counting_routine(A, window, "A")
    C <- counting_routine(C, window, "C")
    G <- counting_routine(G, window, "G")
    U <- counting_routine(U, window, "U")
  }
  counts <- c(A, C, G, U) / length(positions)
  position <- c(rep(seq(1, 10), 4))
  nucleotide <- c(rep("A", 10), rep("C", 10), rep("G", 10), rep("U", 10))
  return(data.frame(counts, position, nucleotide))
}

create_nuc_dist_plot <- function(df, strain, segment, pos) {
  # slice dataset
  df <- df[df$Segment == segment,]
  positions <- df[, pos]

  # load sequence
  sequence <- get_seq(strain, segment)

  # count nuc dist around deletion site
  count_df <- count_nuc_dist(sequence, positions)

  # create a barplot
  ggplot(data=count_df, aes(x=position, y=counts, fill=nucleotide)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_x_continuous(
      breaks=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
      labels=c("5", "4", "3", "2", "1", "-1", "-2", "-3", "-4", "-5")
    )

}

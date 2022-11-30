
counting_routine <- function(l, window, letter, ngs_read_count) {
  count_indices <- unlist(gregexpr(letter, window))
  if (count_indices[1] != -1){
    for (i in count_indices) {
      l[[i]] <- l[[i]] + ngs_read_count
    }
  }
  return(l)
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

create_sampling_data <- function(pos, n_samples, sequence) {
  start <- floor(quantile(pos, probs=seq(0, 1, 1/10))[[2]])
  end <- floor(quantile(pos, probs=seq(0, 1, 1/10))[[10]])
  random_positions <- floor(runif(n_samples, min=start, max=end+1))
  random_counts <- rep(1, n_samples)
  return(count_nuc_dist(sequence, random_positions, random_counts))
}

prepare_nuc_dist_data <- function(df, segment, flattened, strain, sampling) {
  df <- df[df$Segment == segment, ]
  if (nrow(df) == 0) {
    return(df)
  }
  ngs_read_counts <- df[, "NGS_read_count"]
  if (flattened == "flattened") {
    ngs_read_counts[ngs_read_counts != 1] <- 1
  }
 
  # load sequence
  sequence <- get_seq(strain, segment)

  # count nuc dist around deletion site
  if (sampling) {
    n_samples <- nrow(df) * 5
    start_df <- create_sampling_data(df[, "Start"], n_samples, sequence)
    end_df <- create_sampling_data(df[, "End"]-1, n_samples, sequence)
  } else {
    start_df <- count_nuc_dist(sequence, df[, "Start"], ngs_read_counts)
    end_df <- count_nuc_dist(sequence, df[, "End"]-1, ngs_read_counts)
  }
  start_df["location"] <- rep("Start", nrow(start_df))
  end_df["location"] <- rep("End", nrow(end_df))
  count_df <- rbind(start_df, end_df)

  return(count_df)
}

create_nuc_dist_data <- function(df, strain, df2, strain2, segment, flattened){
  # get number of entries or sum of ngs count
  if (flattened == "flattened") {
    ngs_read_count <- nrow(df[df$Segment == segment, ])
  } else {
    ngs_read_count <- sum(df[df$Segment == segment, ]$ngs_read_count)
  }

  count_df1 <- prepare_nuc_dist_data(df, segment, flattened, strain, FALSE)
  # use second data set if availabe
  if (nrow(df2) > 0) {
    count_df1["group"] <- rep("d1", nrow(count_df1))
    count_df2 <- prepare_nuc_dist_data(df2, segment, flattened, strain2, FALSE)
    count_df2["group"] <- rep("d2", nrow(count_df2))
    final_df <- rbind(count_df1, count_df2)
  # create sampling data if no second data set is given
  } else {
    count_df1["group"] <- rep("observed", nrow(count_df1))
    sampling_df <- prepare_nuc_dist_data(df, segment, flattened, strain, TRUE)
    sampling_df["group"] <- rep("expected", nrow(sampling_df))
    final_df <- rbind(count_df1, sampling_df)
  }

  # save as .csv file
  path <- file.path(TEMPPATH, "temp.csv")
  write.csv(final_df, path)
  cat(ngs_read_count, file=file.path(TEMPPATH, "temp.txt"), sep="\n")
}

create_nuc_dist_plot <- function(pos, nuc, segment) {
  # load df and data set length from temp files
  path <- file.path(TEMPPATH, "temp.csv")
  df <- read.csv(path)
  path <- file.path(TEMPPATH, "temp.txt")
  n <- strtoi(readLines(path))
  validate_plotting(df, segment)

  # slice dataset by Start/End and nucleotide
  df <- df[df$location == pos,]
  df <- df[df$nucleotide == nuc,]

  # statistical testing with binom test
  position <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  g1 <- unique(df[c("group")])[[1]][1]
  g2 <- unique(df[c("group")])[[1]][2]
  x_df <- df[df$group == g1,]
  p_df <- df[df$group == g2,]
  p_values <- list()
  for (i in position) {
    x <- x_df[i, "rel_occurrence"] * n
    p <- p_df[i, "rel_occurrence"]
    p_values[[i]] <- binom.test(x, n, p)$p.value
  }
  symbols <- lapply(p_values, get_stat_symbol)

  # max of expected and observed -> is y location of text of stat test
  y_text <- tapply(df$rel_occurrence, df$position, max)
  y_max <- max(0.8, max(y_text)) + 0.05
  # get x coordinates for annotation
  x1 <- ifelse(pos == "Start", 8, 3)
  x2 <- ifelse(pos == "Start", 3, 8)
  # assign x coordinates for grey shadow that marks deletion site
  x_min <- ifelse(pos == "Start", 5.5, 0)
  x_max <- ifelse(pos == "Start", 11, 5.5)
  # create position labels
  labels <- c("5", "4", "3", "2", "1", "-1", "-2", "-3", "-4", "-5")
  if (pos == "End") {
    labels <- rev(labels)
  }

  # create a barplot
  p <- ggplot(
    data=df,
    aes(x=position, y=rel_occurrence, fill=nucleotide, alpha=group)
  ) +
    geom_bar(stat="identity",
      fill=COLOR_MAP[[nuc]],
      color="black",
      position=position_dodge()
    ) +
    ylim(0, y_max) +
    xlab("Position") +
    ylab("Relative occurrence") +
    scale_x_continuous(breaks=position, labels=labels) +
    annotate("text", x=position, y=y_text, label=symbols) +
    geom_rect(xmin=x_min, xmax=x_max, ymin=-1, ymax=2, alpha=0.5, fill="grey")+
    annotate("text", x=x1, y=y_max, label="deleted sequence") +
    annotate("text", x=x2, y=y_max, label="remaining sequence") +
    labs(title=NUC_MAP[[nuc]])

  ggplotly(p)
}


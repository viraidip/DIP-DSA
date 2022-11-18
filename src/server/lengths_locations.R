
create_locations_plot <- function(df, strain, segment, flattened) {
  # slice df by segment, reformat and bind on position and NGS count
  df <- df[df$Segment == segment, ]
  start_df <- df[, c("Start", "NGS_read_count")]
  start_df["Class"] <- "Start"
  colnames(start_df)[colnames(start_df) == "Start"] <- "Position"
  end_df <- df[, c("End", "NGS_read_count")]
  end_df["Class"] <- "End"
  colnames(end_df)[colnames(end_df) == "End"] <- "Position"
  df <- rbind(start_df, end_df)

  # reformat if data should be displayed flattened
  if (flattened == "flattened") {
    df["NGS_read_count"] <- 1
  }

  p <- ggplot(df, aes(x=Position, y=NGS_read_count, fill=Class)) +
    geom_bar(stat="identity") +
    xlab("Nucleotide position on segment") +
    ylab("NGS read count")
  # add info about packaging signal if it exists
  if (packaging_signal_data_exists(strain)) {
    packaging_signal <- load_packaging_signal_data(strain)
    slen <- get_seq_len(strain, segment)
    x <- unlist(packaging_signal[segment])
    y <- layer_scales(p)$y$get_limits()[2]
    color <- c("blue", "blue", "red", "red")
    fill <- c("incorporation signal", "bundling signal")
    p <- p + geom_vline(xintercept=x, color=color, linetype="dotted") +
      geom_rect(aes(xmin=0   , xmax=x[1], ymin=0, ymax=y, fill=fill[1]), alpha=0.3) +
      geom_rect(aes(xmin=x[2], xmax=slen, ymin=0, ymax=y, fill=fill[1]), alpha=0.3) +
      geom_rect(aes(xmin=x[3], xmax=x[1], ymin=0, ymax=y, fill=fill[2]), alpha=0.3) +
      geom_rect(aes(xmin=x[4], xmax=x[2], ymin=0, ymax=y, fill=fill[2]), alpha=0.3)
  }
  ggplotly(p)
}

create_lengths_plot <- function(df, segment, strain, flattened, n_bins) {
  # slice df by segment, reformat and bind on position and NGS count
  df <- df[df$Segment == segment, ]

  seq_len <- get_seq_len(strain, segment)
  df["Length"] <- df["Start"] + (seq_len - df["End"] + 1)
  # multiply each column by NGS count if data is unflattened
  if (flattened != "flattened") {
    df <- data.frame(lapply(df, rep, df$NGS_read_count))
  }

  pl <- ggplot(df, aes(x=Length)) +
    geom_histogram(binwidth=n_bins) +
    xlab("Length of DI candidate") +
    ylab("Number of occurrences")
  # add mean and median to plot
  mean <- mean(df$Length)
  median <- median(df$Length)
  mean_l <- paste("Mean =", format(mean, digits=5))
  median_l <- paste("Median =", format(median, digits=5))
  y <- max(ggplot_build(pl)$data[[1]]$count)
  pl <- pl +
    geom_vline(xintercept=mean, col="#0072B2") +
    annotate("text", x=mean*1.5, y=y, label=mean_l, col="#0072B2") +
    geom_vline(xintercept=median, col="#009E73") +
    annotate("text", x=median*1.5, y=y*0.9, label=median_l, col="#009E73")
  ggplotly(pl)
}


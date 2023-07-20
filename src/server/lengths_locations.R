
format_dataframe_lengths <- function(df, segment, strain, flattened) {
  df <- df[df$Segment == segment, ]
  seq_len <- get_seq_len(strain, segment)
  df["Length"] <- df["Start"] + (seq_len - df["End"] + 1)
  # multiply each column by NGS count if data is unflattened
  if (flattened != "flattened") {
    df <- data.frame(lapply(df, rep, df$NGS_read_count))
  }
  return(df)
}

add_stats <- function(df, pl, class) {
  # select parameters by class
  if (class == "dataset 1") {
    col1 = "#8B2323"
    col2 = "#FF4040"
    y_f = 1.0
  } else {
    col1 = "#27408B"
    col2 = "#4876FF"
    y_f = 0.8
  }
  # calculate stats and add them to plot
  mean <- mean(df$Length)
  median <- median(df$Length)
  mean_l <- paste("Mean =", format(mean, digits=5))
  median_l <- paste("Median =", format(median, digits=5))
  y <- max(ggplot_build(pl)$data[[1]]$count)
  pl <- pl +
    geom_vline(xintercept=mean, col=col1) +
    annotate("text", x=mean, y=y*y_f, label=mean_l, col=col1) +
    geom_vline(xintercept=median, col=col2) +
    annotate("text", x=median, y=y*(y_f-0.1), label=median_l, col=col2)
  return(pl)
}

create_lengths_plot<-function(df,strain,df2,strain2,segment,flattened,n_bins) {
  # slice df by segment, reformat and bind on position and NGS count
  df <- format_dataframe_lengths(df, segment, strain, flattened)
  validate_plotting(df, segment)
  df$Class <- "dataset 1"

  if (nrow(df2) > 0) {
    df2 <- format_dataframe_lengths(df2, segment, strain, flattened)
    if (nrow(df2) > 0) {
      df2$Class <- "dataset 2"
      df <- rbind(df, df2)
    }
  }

  pl <- ggplot(df, aes(x=Length, fill=Class)) +
    geom_histogram(alpha=0.3, binwidth=n_bins, position="identity") +
    xlab("Length of DI candidate") +
    ylab("Number of occurrences") +
    ggtitle(paste("Histogram of DI RNA candidate lengths for segment",
      segment
      )
    ) + 
    theme(plot.title = element_text(size=20))

  # add mean and median to plot
  pl <- add_stats(df[df$Class == "dataset 1", ], pl, "dataset 1")
  if (nrow(df2) > 0) {
    pl <- add_stats(df[df$Class == "dataset 2", ], pl, "dataset 2")
  }
  ggplotly(pl)
}

###############################################################################

format_dataframe_locations <- function(df, segment, flattened) {
  # slice df by segment, reformat and bind on position and NGS count
  df <- df[df$Segment == segment, ]
  if (nrow(df) == 0) {
    return(data.frame())
  }

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

  df2 <- ddply(df,
    c("Position", "Class"),
    summarise,
    NGS_read_count=sum(NGS_read_count)
  )

  return(df2)
}

add_packaging_signal <- function(p, strain, segment) {
  packaging_signal <- load_packaging_signal_data(strain)
  slen <- get_seq_len(strain, segment)
  x <- unlist(packaging_signal[segment])
  y <- layer_scales(p)$y$get_limits()[2]
  color <- c("blue", "blue", "red", "red")
  p <- p + geom_vline(xintercept=x, color=color, linetype="dotted") #+

    '
    geom_rect(
      aes(xmin=0, xmax=x[1], ymin=0, ymax=y, fill="incorporation signal"),
      alpha=0.3
    ) +
    geom_rect(
      aes(xmin=x[2], xmax=slen, ymin=0, ymax=y, fill="incorporation signal"),
      alpha=0.3
    ) +
    geom_rect(
      aes(xmin=x[3], xmax=x[1], ymin=0, ymax=y, fill="bundling signal"),
      alpha=0.3
    ) +
    geom_rect(
      aes(xmin=x[4], xmax=x[2], ymin=0, ymax=y, fill="bundling signal"),
      alpha=0.3
    )
    '
  return (p)
}

create_locations_plot <- function(df, df2, strain, segment, flattened) {
  df <- format_dataframe_locations(df, segment, flattened)
  validate_plotting(df, segment)

  if (nrow(df2) > 0) {
    df2 <- format_dataframe_locations(df2, segment, flattened)
    if (nrow(df2) != 0) {
      df2$Class <- paste(df2$Class, "dataset 2")
      df$Class <- paste(df$Class, "dataset 1")
      df <- rbind(df, df2)
    }
  }

  p <- ggplot(df, aes(x=Position, y=NGS_read_count, fill=Class)) +
    geom_bar(stat="identity", position="dodge", width=1) +
    xlim(0, get_seq_len(strain, segment)) +
    xlab("Nucleotide position on segment") +
    ylab("NGS read count") +
    ggtitle(paste("Start and end position of deletion sites for segment",
      segment
      )
    ) + 
    theme(plot.title = element_text(size=20))

  # add info about packaging signal if it exists
  if (packaging_signal_data_exists(strain)) {
    p <- add_packaging_signal(p, strain, segment)
  }

  ggplotly(p)
}

create_start_end_connection_plot <- function(df, df2, d1, d2, strain, segment, cutoff) {
  df <- df[df$Segment == segment, ]
  validate_plotting(df, segment)
  df["y"] <- 0
  df["yend"] <- 1
  ylab <- d1

  if (nrow(df2) > 0) {
    df2 <- df2[df2$Segment == segment, ]

    if (nrow(df2) > 0) {
      df2["y"] <- 1.2
      df2["yend"] <- 2.2
      df <- rbind(df, df2)
      ylab <- paste(ylab, d2, sep="\t")
    }    
  }

  df <- df[df$NGS_read_count >= cutoff, ]
  max <- max(df$End)

  p <- ggplot() + geom_segment(df,
      mapping=aes(x=Start, y=y, xend=End, yend=yend, col=NGS_read_count)
    ) +
    scale_color_gradient() +
    xlab("Nucleotide position") + 
    ylab(ylab) + 
    geom_rect(aes(xmin=0, xmax=max, ymin=-0.1, ymax=0), alpha=0.9) +
    geom_rect(aes(xmin=0, xmax=max, ymin=1, ymax=1.1), alpha=0.9) +
    annotate(geom="text", x=round(max/2), y=-0.05, label="Start", col="white") +
    annotate(geom="text", x=round(max/2), y=1.05, label="End", col="white")
    ggtitle(paste("Connection of start and end positions for segment",
      segment
      )
    ) + 
    theme(plot.title = element_text(size=20))
  if (nrow(df2) > 0) {
    p <- p + geom_rect(aes(xmin=0, xmax=max, ymin=1.1, ymax=1.2), alpha=0.9) +
    geom_rect(aes(xmin=0, xmax=max, ymin=2.2, ymax=2.3), alpha=0.9) +
    geom_hline(yintercept=1.1, col="red") +
    annotate(geom="text", x=round(max/2), y=1.15, label="Start", col="white") +
    annotate(geom="text", x=round(max/2), y=2.25, label="End", col="white")
  }

  # add info about packaging signal if it exists
  if (packaging_signal_data_exists(strain)) {
    p <- add_packaging_signal(p, strain, segment)
  }

  ggplotly(p)
}

###############################################################################

format_dataframe_end_3_5 <- function(df, strain, segment) {
  df <- df[df$Segment == segment, ]
  validate_plotting(df, segment)
  seq_len <- get_seq_len(strain, segment)
  df["Length_3"] <- df["Start"]
  df["Length_5"] <- seq_len - df["End"]
  return(df)
}

create_end_3_5_plot <- function(df, strain, segment) {
  df <- format_dataframe_end_3_5(df, strain, segment)

  p <- ggplot(df, aes(x=Length_3, y=Length_5)) +
    geom_point() +
    xlim(0, 700) +
    ylim(0, 700) +
    xlab("3' length") +
    ylab("5' length") +
    geom_segment(aes(x=0, xend=900, y=0, yend=900), color="wheat4") +
    
    ggtitle(paste("Relation of 3' and 5' sequence lengths for segment",
      segment
      )
    ) + 
    theme(plot.title = element_text(size=20))
  # add info about packaging signal if it exists
  if (packaging_signal_data_exists(strain)) {
    packaging_signal <- load_packaging_signal_data(strain)
    seq_len <- get_seq_len(strain, segment)
    x <- unlist(packaging_signal[segment])
    l <- c("incorporation signal", "bundling signal")

    xv <- c(x[1], x[3])
    xh <- c(seq_len-x[2], seq_len-x[4])
    color <- c("blue", "red")

    p <- p + geom_vline(xintercept=xv, color=color, linetype="dotted") +
     geom_hline(yintercept=xh, color=color, linetype="dotted")
  }
  ggplotly(p)
}


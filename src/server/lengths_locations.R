
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
  if (class == "1") {
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
  mean_l <- paste("Mean", class,"=", format(mean, digits=5), sep="")
  median_l <- paste("Median", class, "=", format(median, digits=5), sep="")
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
  df$Class <- "1"

  if (nrow(df2) > 0) {
    df2 <- format_dataframe_lengths(df2, segment, strain, flattened)
    df2$Class <- "2"
    df <- rbind(df, df2)
  }

  pl <- ggplot(df, aes(x=Length, fill=Class)) +
    geom_histogram(alpha=0.3, binwidth=n_bins, position="identity") +
    xlab("Length of DI candidate") +
    ylab("Number of occurrences")
  # add mean and median to plot
  pl <- add_stats(df[df$Class == "1", ], pl, "1")
  if (nrow(df2) > 0) {
    pl <- add_stats(df[df$Class == "2", ], pl, "2")
  }
  ggplotly(pl)
}

###############################################################################

format_dataframe_locations <- function(df, segment, flattened) {
  # slice df by segment, reformat and bind on position and NGS count
  df <- df[df$Segment == segment, ]
  validate_plotting(df, segment)

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
  return(df)
}

create_locations_plot <- function(df, df2, strain, segment, flattened) {
  df <- format_dataframe_locations(df, segment, flattened)
  if (nrow(df2) > 0) {
    df2 <- format_dataframe_locations(df2, segment, flattened)
    df2$Class <- paste(df2$Class, "2")
    df <- rbind(df, df2)
  }

  p <- ggplot(df, aes(x=Position, y=NGS_read_count, fill=Class)) +
    geom_bar(stat="identity", width=1) +
    xlim(0, get_seq_len(strain, segment)) +
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
      geom_rect(aes(xmin=0   ,xmax=x[1],ymin=0,ymax=y,fill=fill[1]),alpha=0.3)+
      geom_rect(aes(xmin=x[2],xmax=slen,ymin=0,ymax=y,fill=fill[1]),alpha=0.3)+
      geom_rect(aes(xmin=x[3],xmax=x[1],ymin=0,ymax=y,fill=fill[2]),alpha=0.3)+
      geom_rect(aes(xmin=x[4],xmax=x[2],ymin=0,ymax=y,fill=fill[2]),alpha=0.3)
  }
  ggplotly(p)
}

###############################################################################

format_dataframe_end_3_5 <- function(df, strain, segment) {
  df <- df[df$Segment == segment, ]
  validate_plotting(df, segment)
  seq_len <- get_seq_len(strain, segment)
  df["End_5"] <- seq_len - df["End"]
  return(df)
}

create_end_3_5_plot <- function(df, strain, segment) {
  df <- format_dataframe_end_3_5(df, strain, segment)

  p <- ggplot(df, aes(x=Start, y=End_5)) +
    geom_point() +
    xlim(0, 700) +
    ylim(0, 700) +
    xlab("5' length") +
    ylab("3' length")
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


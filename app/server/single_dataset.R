#### dataset info
plot_ngs_distribution <- function(df1, d1) {
  p <- plot_ly(data = df1, y = ~NGS_read_count, type = "box", name=d1) %>%
    layout(
      title = "Boxplot of NGS read count (Log Scale)",
      yaxis = list(type = "log")
    )

  p
}


#### segment distribution
format_dataframe_segment_distribution <- function(df, flattened) {
  #TODO: add flattened/unlfattened
  segment_counts <- table(df$Segment)
  segment_counts <- segment_counts / sum(segment_counts) * 100
  count_df <- as.data.frame(segment_counts)
  return (count_df)
}

plot_segment_distribution <- function(df, flattened, RCS) {
  df <- apply_cutoff(df, RCS)
  plot_df <- format_dataframe_segment_distribution(df, flattened)
  pl <- ggplot(plot_df, aes(x = Freq, y = factor("All Segments"), fill = Var1)) +
    geom_bar(stat = "identity", color = "white", position = "stack") +
    theme_minimal() +
    labs(x = "Segment of DVG [%]", y = NULL) +
    ggtitle("Segment where DVGs are originating from") + 
    theme(plot.title = element_text(size=20))

  ggplotly(pl)
}


#### deletion shift
format_dataframe_deletion_shift <- function(df, flattened) {
  #TODO: add flattened/unfalltened
  df$del_length <- (df$End-1) - df$Start
  df$shift <- df$del_length %% 3
  shift_counts <- table(df$shift)
  shift_counts <- shift_counts / sum(shift_counts) * 100
  count_df <- as.data.frame(shift_counts)
  return (count_df)
}

plot_deletion_shift <- function(df, flattened, RCS) {
  df <- apply_cutoff(df, RCS)
  plot_df <- format_dataframe_deletion_shift(df, flattened)

  pl <- ggplot(plot_df, aes(x = "", y = Freq, fill = factor(Var1))) +
    geom_bar(stat = "identity", color = "white") +
    ggtitle("Distribution of the deletion shift") +
    theme(plot.title = element_text(size=20))
    labs(fill = "Shift")

  ggplotly(pl)
}


#### lengths and locations
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

plot_lengths<-function(df,strain,segment,flattened,n_bins, RCS) {
  # slice df by segment, reformat and bind on position and NGS count
  df <- apply_cutoff(df, RCS)
  df <- format_dataframe_lengths(df, segment, strain, flattened)
  validate_plotting(df, segment)
  df$Class <- "dataset 1"

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
  ggplotly(pl)
}

format_dataframe_locations <- function(df, segment, flattened, strain) {
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

  df <- ddply(df,
    c("Position", "Class"),
    summarise,
    NGS_read_count=sum(NGS_read_count)
  )

  get_nucleotide_at_position <- function(position, seq) {
    as.character(subseq(seq, start=position, end=position))
  }
  sequence <- get_seq(strain, segment)
  df$Nucleotide <- sapply(df$Position, get_nucleotide_at_position, seq=sequence)

  return(df)
}

add_packaging_signal <- function(p, strain, segment) {
  packaging_signal <- load_packaging_signal_data(strain)
  slen <- get_seq_len(strain, segment)
  x <- unlist(packaging_signal[segment])
  y <- layer_scales(p)$y$get_limits()[2]
  color <- c("blue", "blue", "red", "red")
  p <- p + geom_vline(xintercept=x, color=color, linetype="dotted") 
  return (p)
}

plot_locations <- function(df, strain, segment, flattened, RCS) {
  df <- apply_cutoff(df, RCS)
  df <- format_dataframe_locations(df, segment, flattened, strain)
  validate_plotting(df, segment)

  p <- ggplot(df, aes(x=Position, y=NGS_read_count, fill=Nucleotide)) +
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

plot_start_end_connection <- function(df, d1, strain, segment, RCS) {
  df <- apply_cutoff(df, RCS)
  df <- df[df$Segment == segment, ]
  validate_plotting(df, segment)
  df["y"] <- 0
  df["yend"] <- 1
  ylab <- d1

  max <- get_seq_len(strain, segment) 

  p <- ggplot() + geom_segment(df,
      mapping=aes(x=Start, y=y, xend=End, yend=yend, col=NGS_read_count)
    ) +
    scale_color_gradient() +
    xlab("Nucleotide position") + 
    ylab(ylab) + 
    geom_rect(aes(xmin=0, xmax=max, ymin=-0.1, ymax=0), alpha=0.9) +
    geom_rect(aes(xmin=0, xmax=max, ymin=1, ymax=1.1), alpha=0.9) +
    annotate(geom="text", x=round(max/2), y=-0.05, label="Start", col="white") +
    annotate(geom="text", x=round(max/2), y=1.05, label="End", col="white") +
    ggtitle(paste("Connection of start and end positions for segment",
      segment
      )
    ) + 
    theme(plot.title = element_text(size=20))

  ggplotly(p)
}

format_dataframe_end_3_5 <- function(df, strain, segment) {
  df <- df[df$Segment == segment, ]
  validate_plotting(df, segment)
  seq_len <- get_seq_len(strain, segment)
  df["Length_3"] <- df["Start"]
  df["Length_5"] <- seq_len - df["End"]
  return(df)
}

plot_end_3_5_plot <- function(df, strain, segment, RCS) {
  df <- apply_cutoff(df, RCS)
  df <- format_dataframe_end_3_5(df, strain, segment)

  p <- ggplot(df, aes(x=Length_3, y=Length_5)) +
    geom_point() +
    geom_abline(a=0, b=1, col="green") +
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

#### nucleotide enrichment
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

prepare_nuc_dist_data <- function(df, segment, flattened, strain) {
  df <- df[df$Segment == segment, ]
  if (nrow(df) == 0) {
    return(df)
  }
  
  if (flattened == "flattened") {
    df$NGS_read_count <- 1
  }
  ngs_read_counts <- df$NGS_read_count

  # load sequence
  sequence <- get_seq(strain, segment)

  # count nuc dist around deletion site
  start_df <- count_nuc_dist(sequence, df[, "Start"], ngs_read_counts)
  end_df <- count_nuc_dist(sequence, df[, "End"]-1, ngs_read_counts)
  start_df["location"] <- rep("Start", nrow(start_df))
  end_df["location"] <- rep("End", nrow(end_df))
  count_df <- rbind(start_df, end_df)
  return(count_df)
}

plot_nuc_dist <- function(pos, nuc, segment, datasetname, RCS, strain, flattened) {
  # load df and dataset length from temp files
  df <- read.csv(file.path(DATASETSPATH, strain, paste(datasetname, ".csv", sep="")), na.strings=c("NaN"))
  exp_df <- read.csv(file.path(DATASETSPATH, strain, paste(datasetname, ".tsv", sep="")), na.strings=c("NaN"), sep="\t")

  df <- apply_cutoff(df, RCS)

  # get number of entries or sum of ngs count
  if (flattened == "flattened") {
    ngs_read_count <- nrow(df)
  } else {
    ngs_read_count <- sum(df$NGS_read_count)
  }

  count_df <- prepare_nuc_dist_data(df, segment, flattened, strain)
  count_df["group"] <- rep("observed", nrow(count_df))
  n <- nrow(count_df)
  sampling_df <- prepare_nuc_dist_data(exp_df, segment, flattened, strain)
  sampling_df["group"] <- rep("expected", nrow(sampling_df))
  final_df <- rbind(count_df, sampling_df)

  validate_plotting(final_df, segment)

  # slice dataset by Start/End and nucleotide
  final_df <- final_df[final_df$location == pos,]
  final_df <- final_df[final_df$nucleotide == nuc,]

  # statistical testing with binom test
  position <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  g1 <- unique(final_df[c("group")])[[1]][1]
  g2 <- unique(final_df[c("group")])[[1]][2]
  x_df <- final_df[final_df$group == g1,]
  p_df <- final_df[final_df$group == g2,]
  p_values <- list()
  for (i in position) {
    x <- as.integer(x_df[i, "rel_occurrence"] * n)
    p <- p_df[i, "rel_occurrence"]
    p_values[[i]] <- binom.test(x, n, p)$p.value
  }
  symbols <- lapply(p_values, get_stat_symbol)

  # max of expected and observed -> is y location of text of stat test
  y_text <- tapply(final_df$rel_occurrence, final_df$position, max)
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
    data=final_df,
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
    ggtitle(paste("Nucleotide distribution of",
      NUC_MAP[[nuc]],
      "at",
      pos,
      "for segment",
      segment
      )
    ) + 
    theme(plot.title = element_text(size=20)) +
    scale_x_continuous(breaks=position, labels=labels) +
    annotate("text", x=position, y=y_text, label=symbols) +
    geom_rect(xmin=x_min, xmax=x_max, ymin=-1, ymax=2, alpha=0.5, fill="grey")+
    annotate("text", x=x1, y=y_max, label="deleted sequence") +
    annotate("text", x=x2, y=y_max, label="remaining sequence")

  ggplotly(p)
}

#### direct repeats
direct_repeats_counting_routine <- function(row, sequence) {
  start <- as.integer(row["Start"])
  end <- as.integer(row["End"])
  i_s_1 <- start-5
  i_s_2 <- start
  i_e_1 <- end-5-1
  i_e_2 <- end-1
  start_window <- sequence[i_s_1:i_s_2]
  end_window <- sequence[i_e_1:i_e_2]
  counter <- 0
  for (i in length(start_window):1) {
    if (start_window[i] == end_window[i]) {
      counter <- counter + 1
    } else {
      return(counter)
    }
  }

  return(counter)
}


prepare_plot_data <- function(df, label) {
  # get data set by label
  df <- df[df$group == label, ]

  # calculate direct repeat lengths as ratio
  table <- table(df$direct_repeats)
  table <- table/sum(table)

  df <- data.frame(table, rep(label, length(table)))
  colnames(df) <- c("length", "freq", "group")
  df$length <- as.numeric(as.character(df$length))

  return(df)
}

plot_direct_repeats <- function(strain, datasetname, segment, RCS, flattened) {
  df <- read.csv(file.path(DATASETSPATH, strain, paste(datasetname, ".csv", sep="")), na.strings=c("NaN"))
  exp_df <- read.csv(file.path(DATASETSPATH, strain, paste(datasetname, ".tsv", sep="")), na.strings=c("NaN"), sep="\t")

  df <- apply_cutoff(df, RCS)
  df <- df[df$Segment == segment,]
  exp_df <- exp_df[exp_df$Segment == segment,]
  df <- subset(df, select=-c(Segment))
  exp_df <- subset(exp_df, select=-c(Segment))

  # include NGS count or not
  if (flattened == "flattened") {
    df$NGS_read_count[df$NGS_read_count != 1] <- 1
    r_df <- df
  } else {
    r_df <- data.frame(Start=integer(),End=integer(),NGS_read_count=integer())
    for (i in 1:nrow(df)) {
      r_df <- rbind(r_df, df[rep(i, df[i,3]),])
    }
  }

  if (nrow(r_df) == 0) {
    return()
  }
  
  seq <- get_seq(strain, segment)
  n_samples <- nrow(r_df)

  r_df["group"] <- rep("observed", n_samples)
  exp_df["group"] <- rep("expected", nrow(exp_df))

  final_df <- rbind(r_df, exp_df)
  final_df["direct_repeats"] <- apply(final_df,
    1,
    direct_repeats_counting_routine,
    seq
  )
  
  validate_plotting(final_df, segment)

  g1 <- unique(final_df[c("group")])[[1]][1]
  g2 <- unique(final_df[c("group")])[[1]][2]
  df_1 <- prepare_plot_data(final_df, g1)
  df_2 <- prepare_plot_data(final_df, g2)
  n_1 <- nrow(df_1)
  n_2 <- nrow(df_2)

  # fill NA values with 0.0 to have a good representation in the final plot
  max_length <- max(max(df_1$length), max(df_2$length))
  for (i in 0:max_length) {
    if (!any(i==df_1[,1])) {
      df_1 <- rbind(df_1, list(i, 0.0, g1))
    }
    if (!any(i==df_2[,1])) {
      df_2 <- rbind(df_2, list(i, 0.0, g2))
    }
  }

  plot_df <- merge(df_1, df_2, all=TRUE)
  plot_df$freq <- as.numeric(plot_df$freq)

  # statistical testing with Wilcoxon/Mann-Witney test

  data_1 <- df[df$group == g1, ]$direct_repeats
  data_2 <- df[df$group == g2, ]$direct_repeats
  

  if (length(data_1) == 0) {
    symbol <- ""
  } else {
    res <- wilcox.test(data_1, data_2)
    symbol <- get_stat_symbol(res$p.value)
  }
    
  t1 <- "Frequency of different direct repeat lengths for segment "
  t2 <- paste(segment, " (n=", n_samples, ") ", symbol, sep="")
  title <- paste(t1, t2, sep="")

  # create a barplot
  p <- ggplot(data=plot_df, aes(x=length, y=freq, fill=group)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ylim(0, 1.0) +
    xlab("Length of direct repeat") +
    ylab("Relative occurrence") +
    ggtitle(title) +
    theme(plot.title = element_text(size=20))
  ggplotly(p)
}


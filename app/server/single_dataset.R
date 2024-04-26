plot_ngs_distribution <- function(strain, datasetname, RSC, prg) {
  df <- load_single_dataset(file.path(strain,paste(datasetname,".csv",sep="")))
  df <- apply_cutoff(df, RSC)
  validate_df(df)
  label <- paste(datasetname, " (n=", nrow(df), ")", sep="")
  df$name <- datasetname

  pl <- ggplot(df, aes(x=name, y=NGS_read_count, fill=name)) +
    geom_boxplot() +
    scale_y_log10(limits=c(1, NA)) +
    labs(x=label, y="NGS count (log scale)") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  prg$set(0.1, "NGS count plot")
  ggplotly(pl)
}


plot_frame_shift <- function(strain, datasetname, flattened, RSC, prg) {
  df <- load_single_dataset(file.path(strain,paste(datasetname,".csv",sep="")))
  df <- apply_cutoff(df, RSC)
  validate_df(df)

  if (flattened == "flattened") {
    df["NGS_read_count"] <- 1
  }

  df$del_length <- (df$End-1) - df$Start
  df$Shift <- df$del_length %% 3

  plot_df <- df %>%
    group_by(Shift) %>%
    summarise(counts=sum(NGS_read_count)) %>%
    mutate(Freq=counts / sum(counts) * 100)

  plot_df <- plot_df %>%
    mutate(Shift=case_when(
      Shift == 0 ~ "in-frame",
      Shift == 1 ~ "shift +1",
      Shift == 2 ~ "shift -1",
      TRUE ~ as.character(Shift)
    ))

  expected <- c(1/3, 1/3, 1/3)
  r <- chisq.test(plot_df[c("Freq")], p=expected)
  label <- paste(datasetname, get_stat_symbol(r$p.value))

  pl <- ggplot(plot_df, aes(x="", y=Freq, fill=factor(Shift))) +
    geom_bar(stat="identity") +
    labs(x=label, y="Deletion shift [%]", fill="Shifts")
  
  prg$set(0.2, "Frame shift plot")
  ggplotly(pl)
}


plot_segment_distribution <- function(strain, datasetname, flattened, RSC, prg) {
  df <- load_single_dataset(file.path(strain,paste(datasetname,".csv",sep="")))
  df <- apply_cutoff(df, RSC)
  validate_df(df)

  if (flattened == "flattened") {
    df["NGS_read_count"] <- 1
  }

  plot_df <- df %>%
    group_by(Segment) %>%
    summarise(counts=sum(NGS_read_count)) %>%
    mutate(Freq=counts / sum(counts) * 100)
  
  plot_df <- plot_df %>%
    rowwise() %>%
    mutate(seq_len=get_seq_len(strain, Segment))

  plot_df$expected <- plot_df$seq_len / sum(plot_df$seq_len) * 100
  r <- chisq.test(plot_df[c("Freq", "expected")])
  label <- paste(datasetname, "", get_stat_symbol(r$p.value))

  plot_df$name <- datasetname
  pl <- ggplot(plot_df, aes(x=name, y=Freq, fill=Segment)) +
    geom_bar(stat="identity", position="stack") +
    labs(x=label, y="Segment of DelVG [%]") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  prg$set(0.3, "Segment distribution plot")
  ggplotly(pl)
}


plot_lengths<-function(strain, datasetname, segment,flattened,n_bins, RSC, prg) {
  df <- load_single_dataset(file.path(strain,paste(datasetname,".csv",sep="")))
  df <- apply_cutoff(df, RSC)
  df <- format_dataframe_lengths(df, segment, strain, flattened)
  validate_plotting(df, segment)
  df$Class <- datasetname

  pl <- ggplot(df, aes(x=Length, fill=Class)) +
    geom_histogram(bins=n_bins, position="identity") +
    labs(x="Length of DelVG", y="Number of occurrences", fill="Dataset")
  # add mean and median to plot
  pl <- add_stats_lengths(df, pl)

  prg$set(0.4, "Deletion lengths plot")
  ggplotly(pl)
}


plot_locations <- function(strain, datasetname, segment, flattened, RSC, prg) {
  df <- load_single_dataset(file.path(strain,paste(datasetname,".csv",sep="")))
  df <- apply_cutoff(df, RSC)
  df <- format_dataframe_locations(df, segment, flattened, strain)
  validate_plotting(df, segment)

  if (flattened == "flattened") {
    y_lab <- "Number of occurrences"
  } else {
    y_lab <- "NGS read count"
  }

  p <- ggplot(df, aes(x=Position, y=NGS_read_count, fill=Nucleotide)) +
    geom_bar(stat="identity", position="dodge", width=1) +
    xlim(0, get_seq_len(strain, segment)) +
    labs(x="Nucleotide position on segment", y=y_lab)

  # add info about packaging signal if it exists
  if (packaging_signal_data_exists(strain)) {
    p <- add_packaging_signal(p, strain, segment)
  }

  prg$set(0.5, "Deletion location plot")
  ggplotly(p)
}


plot_end_3_5 <- function(strain, datasetname, segment, RSC, prg) {
  df <- load_single_dataset(file.path(strain,paste(datasetname,".csv",sep="")))
  df <- apply_cutoff(df, RSC)
  
  df <- df[df$Segment == segment, ]
  validate_plotting(df, segment)
  seq_len <- get_seq_len(strain, segment)
  df["Length_3"] <- df["Start"]
  df["Length_5"] <- seq_len - df["End"]

  p <- ggplot(df, aes(x=Length_3, y=Length_5)) +
    geom_point() +
    geom_abline(intercept=0, slope=1, col="#00BA38") +
    xlim(0, 700) +
    ylim(0, 700) +
    labs(x="3' length", y="5' length")

  # add info about packaging signal if it exists
  if (packaging_signal_data_exists(strain)) {
    packaging_signal <- load_packaging_signal_data(strain)
    seq_len <- get_seq_len(strain, segment)
    x <- unlist(packaging_signal[segment])
    l <- c("incorporation signal", "bundling signal")

    xv <- c(x[1], x[3])
    xh <- c(seq_len-x[2], seq_len-x[4])
    color <- c("#619CFF", "#F8766D")

    p <- p + geom_vline(xintercept=xv, color=color, linetype="dotted") +
     geom_hline(yintercept=xh, color=color, linetype="dotted")
  }

  prg$set(0.6, "3' 5' plot")
  ggplotly(p)
}


plot_start_end_mapping <- function(strain, datasetname, segment, RSC, prg) {
  df <- load_single_dataset(file.path(strain,paste(datasetname,".csv",sep="")))
  df <- apply_cutoff(df, RSC)
  df <- df[df$Segment == segment, ]
  validate_plotting(df, segment)
  df["y"] <- 0
  df["yend"] <- 1

  max <- get_seq_len(strain, segment) 

  p <- ggplot() +
    geom_segment(
      df,
      mapping=aes(x=Start, y=y, xend=End, yend=yend, color=NGS_read_count)
    ) +
    scale_color_gradient() +
    labs(x="Nucleotide position", y=datasetname, color="NGS read count") +
    geom_rect(aes(xmin=0, xmax=max, ymin=-0.1, ymax=0), alpha=0.9) +
    geom_rect(aes(xmin=0, xmax=max, ymin=1, ymax=1.1), alpha=0.9) +
    annotate(geom="text", x=round(max/2), y=-0.05, label="Start", col="white")+
    annotate(geom="text", x=round(max/2), y=1.05, label="End", col="white")

  prg$set(0.7, "Start and end plot")
  ggplotly(p)
}


plot_direct_repeats <- function(strain, datasetname, segment, RSC, flattened, prg) {
  df <- load_single_dataset(
    file.path(strain, paste(datasetname, ".csv", sep=""))
  )
  exp_df <- load_single_dataset(
    file.path(strain, paste(datasetname, ".tsv", sep="")),
    sep="\t"
  )

  df <- apply_cutoff(df, RSC)
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

  validate_plotting(r_df, segment)
  
  seq <- get_seq(strain, segment)
  n_samples <- nrow(r_df)
  r_df$direct_repeats <- apply(r_df,1,direct_repeats_counting_routine,seq)
  df_1 <- prepare_direct_repeat_plot_data(r_df, "observed")
  
  # only calculate direct repeats if expected data is available
  if (nrow(exp_df) != 0) {
    exp_df$direct_repeats <- apply(exp_df,1,direct_repeats_counting_routine,seq)
    df_2 <- prepare_direct_repeat_plot_data(exp_df, "expected")
    do_testing <- TRUE
  } else {
    df_2 <- data.frame(matrix(ncol=length(names(df_1)), nrow=0))
    names(df_2) <- names(df_1)
    do_testing <- FALSE
  }


  # fill NA values with 0.0 to have a good representation in the final plot
  max_length <- max(max(df_1$length), max(df_2$length))
  for (i in 0:max_length) {
    if (!any(i==df_1[,1])) {
      df_1 <- rbind(df_1, list(i, 0, "observed", 0.0))
    }
    if (!any(i==df_2[,1])) {
      df_2 <- rbind(df_2, list(i, 0, "expected", 0.0))
    }
  }

  plot_df <- merge(df_1, df_2, all=TRUE)
  plot_df$freq <- as.numeric(plot_df$freq)

  # statistical testing with Chi-squared test
  if (do_testing) {
    matrix <- matrix(c(plot_df[plot_df$group == "observed", "counts"],
                    plot_df[plot_df$group == "expected", "counts"]), ncol=2)
    data_1 <- r_df$direct_repeats
    data_2 <- exp_df$direct_repeats
    if (length(data_1) == 0) {
      s <- ""
    } else {
      r <- chisq.test(matrix)
      p <- r$p.value
      s <- get_stat_symbol(p)
    }
  } else {
    s <- ""
  }

  text <- paste("n=", n_samples, ", ", s, sep="")
  # create a barplot
  p <- ggplot(data=plot_df, aes(x=length, y=freq, fill=group)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ylim(0, 1.0) +
    labs(x="Length of direct repeat", y="Relative occurrence",fill="Dataset") +
    annotate("text", x=3, y=0.9, label=text) +
    theme(plot.title = element_text(size=20))
  
  prg$set(0.8, "Direct repeats plot")
  ggplotly(p)
}


plot_nucleotide_enrichment <- function(strain,
                                       datasetname,
                                       pos,
                                       nuc,
                                       segment,
                                       RSC,
                                       flattened,
                                       prg) {  
  df <- load_single_dataset(
    file.path(strain, paste(datasetname, ".csv", sep=""))
  )
  exp_df <- load_single_dataset(
    file.path(strain, paste(datasetname, ".tsv", sep="")),
    sep="\t"
  )

  df <- apply_cutoff(df, RSC)
  n <- sum(df$Segment == segment)
  count_df <- prepare_nuc_enr_data(df, segment, flattened, strain, pos, nuc)
  count_df["group"] <- rep("observed", nrow(count_df))
  
  # only calculate nucleotide enrichment if expected data is available
  sampling_df <- prepare_nuc_enr_data(exp_df,segment,flattened,strain,pos,nuc)
  if (nrow(sampling_df) != 0) {
    sampling_df["group"] <- rep("expected", nrow(sampling_df))
    do_testing <- TRUE
  } else {
    do_testing <- FALSE
  }
  
  final_df <- rbind(count_df, sampling_df)
  validate_plotting(final_df, segment)

  # statistical testing with binom test
  # only performed when sampling data exists for this segment
  position <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  if (do_testing) {
    p_values <- list()
    for (i in position) {
      x <- as.integer(count_df[i, "rel_occurrence"] * n)
      p <- sampling_df[i, "rel_occurrence"]
      p_values[[i]] <- binom.test(x, n, p)$p.value
    symbols <- lapply(p_values, get_stat_symbol)
    symbols <- gsub("ns.", "", symbols)
    }
  } else {
    symbols <- rep("", 10)
  }
  
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
      color=COLOR_MAP[[nuc]],
      position=position_dodge()
    ) +
    ylim(0, y_max) +
    labs(x="Position", y="Relative occurrence", alpha="Dataset") +
    theme(plot.title=element_text(size=20)) +
    scale_x_continuous(breaks=position, labels=labels) +
    annotate("text", x=position, y=y_text+0.02, label=symbols, size=6) +
    geom_rect(xmin=x_min, xmax=x_max, ymin=-1, ymax=2, alpha=0.5, fill="grey")+
    annotate("text", x=x1, y=y_max, label="deleted sequence") +
    annotate("text", x=x2, y=y_max, label="remaining sequence")

  if (pos == "Start") {
    prg$set(0.9, "Nucleotide enrichment plot")
  } else {
    prg$close()
  }
  ggplotly(p)
}


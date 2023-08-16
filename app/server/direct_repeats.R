
direct_repeats_counting_routine <- function(row, sequence) {
  start <- as.integer(row["Start"])
  end <- as.integer(row["End"])
  i_s_1 <- start-15
  i_s_2 <- start
  i_e_1 <- end-15-1
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

generate_positions <- function(pos, n_samples) {
  start <- floor(quantile(pos, probs=seq(0, 1, 1/10))[[2]])
  end <- floor(quantile(pos, probs=seq(0, 1, 1/10))[[10]])
  random_positions <- floor(runif(n_samples, min=start, max=end+1))
  return(random_positions)
}

create_direct_repeat_sampling_data <- function(df, n_samples, sequence) {
  Start <- generate_positions(df[, "Start"], n_samples)
  End <- generate_positions(df[, "End"], n_samples)
  NGS_read_count <- rep(1, n_samples)

  return(data.frame(Start, End, NGS_read_count))
}

prepare_data <- function(df, strain, segment, flattened) {
  # load observed data
  df <- df[df$Segment == segment,]
  df <- subset(df, select=-c(Segment))

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

  return(r_df)
}

create_direct_repeats_data<-function(df,strain,df2,strain2,segment,flattened) {
  df <- prepare_data(df, strain, segment, flattened)
  if (nrow(df) == 0) {
    path <- file.path(TEMPPATH, "direct_repeats_temp.csv")
    write.csv(df, path)
    path <- file.path(TEMPPATH, "direct_repeats_temp.txt")
    cat(0, file=path, sep="\n")
    return()
  }
  s <- get_seq(strain, segment)

  # check if a second dataset is given to compare
  if (nrow(df2) > 0) {
    df2 <- prepare_data(df2, strain2, segment, flattened)
    # catch error if second dataset has no entries on selected strain
    if (nrow(df2) != 0) {
      df["direct_repeats"] <- apply(df, 1, direct_repeats_counting_routine, s)
      df["group"] <- rep("dataset 1", nrow(df))
      s2 <- get_seq(strain2, segment)
      df2["direct_repeats"] <- apply(df2, 1, direct_repeats_counting_routine, s2)
      df2["group"] <- rep("dataset 2", nrow(df2))

      final_df <- rbind(df, df2)
    } else {
      df["group"] <- rep("observed", nrow(df))
      n_samples <- nrow(df) * 5
      sampling_df <- create_direct_repeat_sampling_data(df, n_samples, s)
      sampling_df["group"] <- rep("expected", nrow(sampling_df))

      final_df <- rbind(df, sampling_df)
      final_df["direct_repeats"] <- apply(final_df,
        1,
        direct_repeats_counting_routine,
        s
      )
    }
  } else {
    df["group"] <- rep("observed", nrow(df))
    # create sampling data
    n_samples <- nrow(df) * 5
    sampling_df <- create_direct_repeat_sampling_data(df, n_samples, s)
    sampling_df["group"] <- rep("expected", nrow(sampling_df))

    final_df <- rbind(df, sampling_df)
    final_df["direct_repeats"] <- apply(final_df,
      1,
      direct_repeats_counting_routine,
      s
    )
  }

  # save as .csv file
  path <- file.path(TEMPPATH, "direct_repeats_temp.csv")
  write.csv(final_df, path)
  path <- file.path(TEMPPATH, "direct_repeats_temp.txt")
  cat(sum(df$NGS_read_count), file=path, sep="\n")
}

add_correction <- function(df) {
  df["freq_c"] <- rep(0.0, nrow(df))
  for (i in df$length) {
    orig_value <- df[df$length == i, "freq"]
    if (orig_value != 0) {
      divisor <- orig_value/(i+1)
      df[df$length == i, "freq_c"] <- divisor
      for (j in 0:i-1) {
        df[df$length == j, "freq_c"] <- df[df$length == j, "freq_c"] + divisor
      }
    }
  }

  df$freq <- df$freq_c
  df <- subset(df, select=-c(freq_c))
  return(df)
}

prepare_plot_data <- function(df, label, correction) {
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

create_direct_repeats_plot <- function(correction, segment) {
  # load df and data set length from temp files
  path <- file.path(TEMPPATH, "direct_repeats_temp.csv")
  df <- read.csv(path)
  path <- file.path(TEMPPATH, "direct_repeats_temp.txt")
  n_samples <- strtoi(readLines(path))
  validate_plotting(df, segment)

  g1 <- unique(df[c("group")])[[1]][1]
  g2 <- unique(df[c("group")])[[1]][2]

  df_1 <- prepare_plot_data(df, g1, correction)
  # df_2 is either expected data or second data set if one is selected
  df_2 <- prepare_plot_data(df, g2, correction)
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

  # add correction factor if wanted
  if (correction == "Yes") {
    df_1 <- add_correction(df_1)
    # also add two second set, if not expected values
    if (g2 == "d2") {
      df_2 <- add_correction(df_2)
    }
  }

  plot_df <- merge(df_1, df_2, all=TRUE)
  plot_df$freq <- as.numeric(plot_df$freq)

  # statistical testing with Wilcoxon/Mann-Witney test
  if (correction == "Yes") {
    data_1 <- c()
    data_2 <- c()
    for (i in 0:max_length) {
      n <- plot_df[plot_df$length == i & plot_df$group == g1, "freq"]
      data_1 <- c(data_1, rep(i, round(n*n_1)))
      n <- plot_df[plot_df$length == i & plot_df$group == g2, "freq"]
      data_2 <- c(data_2, rep(i, round(n*n_2)))
    }
  }
  else {
    data_1 <- df[df$group == g1, ]$direct_repeats
    data_2 <- df[df$group == g2, ]$direct_repeats
  }

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


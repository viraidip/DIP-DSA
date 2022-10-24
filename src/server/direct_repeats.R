
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
}

create_direct_repeat_sampling_data <- function(df, n_samples, sequence) {
  Start <- generate_positions(df[, "Start"], n_samples)
  End <- generate_positions(df[, "End"], n_samples)
  NGS_read_count <- rep(1, n_samples)

  data.frame(Start, End, NGS_read_count)
}

create_direct_repeats_data <- function(df, strain, segment, flattened) {
  # load observed data
  df <- df[df$Segment == segment,]
  df <- subset(df, select=-c(Segment))
  if (flattened == "flattened") {
    df$NGS_read_count[df$NGS_read_count != 1] <- 1
  }

  # load sequence
  sequence <- get_seq(strain, segment)

  # count nuc dist around deletion site
  count_df <- df
  count_df["group"] <- rep("observed", nrow(count_df))
 
  # create sampling data
  n_samples <- nrow(df) * 5
  sampling_df <- create_direct_repeat_sampling_data(df, n_samples, sequence)
  sampling_df["group"] <- rep("expected", nrow(sampling_df))

  final_df <- rbind(count_df, sampling_df)
  final_df["direct_repeats"] <- apply(final_df, 1, direct_repeats_counting_routine, sequence)

  # save as .csv file
  path <- file.path(TEMPPATH, "direct_repeats_temp.csv")
  write.csv(final_df, path)
  
  cat(sum(df$NGS_read_count), file=file.path(TEMPPATH, "direct_repeats_temp.txt"), sep="\n")
}

add_correction <- function(df) {
  df["freq_corrected"] <- rep(0.0, nrow(df))
  for (i in df$length) {
    orig_value <- df[df$length == i, "freq"]
    if (orig_value != 0) {
      divided_value <- orig_value/(i+1)
      df[df$length == i, "freq_corrected"] <- divided_value
      for (j in 0:i-1) {
        df[df$length == j, "freq_corrected"] <- df[df$length == j, "freq_corrected"] + divided_value
      }
    }
  }

  df$freq <- df$freq_corrected
  df <- subset(df, select=-c(freq_corrected))
  return(df)
}

create_direct_repeats_plot <- function(correction) {
  # load df and data set length from temp files
  path <- file.path(TEMPPATH, "direct_repeats_temp.csv")
  df <- read.csv(path)
  path <- file.path(TEMPPATH, "direct_repeats_temp.txt")
  n_samples <- strtoi(readLines(path))

  df$direct_repeats <- df$direct_repeats * df$NGS_read_count
  obs_df <- df[df$group == "observed", ]
  exp_df <- df[df$group == "expected", ]

  n_obs <- nrow(obs_df)
  n_exp <- nrow(exp_df)

  obs_table <- table(obs_df$direct_repeats)
  obs_table <- obs_table/sum(obs_table)
  exp_table <- table(exp_df$direct_repeats)
  exp_table <- exp_table/sum(exp_table)

  obs_df <- data.frame(obs_table, rep("observed", length(obs_table)))
  colnames(obs_df) <- c("length", "freq", "group")
  obs_df$length <- as.numeric(as.character(obs_df$length))
  exp_df <- data.frame(exp_table, rep("expected", length(exp_table)))
  colnames(exp_df) <- c("length", "freq", "group")
  exp_df$length <- as.numeric(as.character(exp_df$length))

  max_length <- max(max(obs_df$length), max(exp_df$length))
  for (i in 0:max_length) {
    if (!any(i==obs_df[,1])) {
      obs_df <- rbind(obs_df, list(i, 0.0, "observed"))
    }
    if (!any(i==exp_df[,1])) {
      exp_df <- rbind(exp_df, list(i, 0.0, "expected"))
    }
  }

  if (correction == "Yes") {
    obs_df <- add_correction(obs_df)
  }

  plot_df <- merge(obs_df, exp_df, all=TRUE)
  plot_df$freq <- as.numeric(plot_df$freq)

  # statistical testing with Wilcoxon/Mann-Witney test
  if (correction == "Yes") {
    obs_data <- c()
    exp_data <- c()
    for (i in 0:max_length) {
      n <- plot_df[plot_df$length == i & plot_df$group == "observed", "freq"]
      obs_data <- c(obs_data, rep(i, round(n*n_obs)))
      n <- plot_df[plot_df$length == i & plot_df$group == "expected", "freq"]
      exp_data <- c(exp_data, rep(i, round(n*n_exp)))
    }
  }
  else {
    obs_data <- df[df$group == "observed", ]$direct_repeats
    exp_data <- df[df$group == "expected", ]$direct_repeats
  }
  res <- wilcox.test(obs_data, exp_data)
  symbol <- get_stat_symbol(res$p.value)

  # create a barplot
  p <- ggplot(data=plot_df, aes(x=length, y=freq, fill=group)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ylim(0, 1.0) +
    ggtitle(paste("frequency of different direct repeat lengths (n=", n_samples, ") ", symbol, sep=""))
  ggplotly(p)
}


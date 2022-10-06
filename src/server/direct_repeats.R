

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

  cat(sum(df$NGS_read_counts), file=file.path(TEMPPATH, "direct_repeats_temp.txt"), sep="\n")
}

create_direct_repeats_plot <- function() {
  # load df and data set length from temp files
  path <- file.path(TEMPPATH, "direct_repeats_temp.csv")
  df <- read.csv(path)
  path <- file.path(TEMPPATH, "direct_repeats_temp.txt")
  n <- strtoi(readLines(path))

  df$direct_repeats <- df$direct_repeats * df$NGS_read_count
  obs_df <- df[df$group == "observed", ]
  exp_df <- df[df$group == "expected", ]

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

  for (i in 1:nrow(obs_df)) {
    element <- obs_df[i,1]
    if (!any(element==exp_df[,1])) { #if number is not in exp_df add to exp_df
      exp_df <- rbind(exp_df, c(element, 0.0, "expected"))
    }
  }
  for (i in 1:nrow(exp_df)) {
    element <- exp_df[i,1]
    if (!any(element==obs_df)) { #if number is not in obs_df add to obs_df
      obs_df <- rbind(obs_df, c(element, 0.0, "observed"))
    }
  }

  plot_df <- merge(obs_df, exp_df, all=TRUE)
  plot_df$freq <- as.numeric(plot_df$freq)

  # statistical testing with Wilcoxon/Mann-Witney test
  obs_data <- df[df$group == "observed", ]$direct_repeats
  exp_data <- df[df$group == "expected", ]$direct_repeats
  res <- wilcox.test(obs_data, exp_data)
  symbol <- get_stat_symbol(res$p.value)

  # create a barplot
  p <- ggplot(data=plot_df, aes(x=length, y=freq, alpha=group)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    annotate("text", x=5, y=0.5, label=symbol)
  ggplotly(p)
}


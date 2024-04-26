plot_multiple_ngs_distribution <- function(paths, RSC, prg) {
  validate_selection(paths)
  df <- load_all_datasets(paths)
  df <- apply_cutoff(df, RSC)
  validate_df(df)

  pl <- ggplot(df, aes(x=name, y=NGS_read_count, fill=name)) +
    geom_boxplot() +
    scale_y_log10(limits=c(1, NA)) +
    labs(x="Dataset", y="NGS count (log scale)", fill="Dataset")

  prg$set(0.1, "NGS count plot")
  ggplotly(pl)
}


plot_multiple_deletion_shift <- function(paths, flattened, RSC, prg) {
  validate_selection(paths)
  df <- load_all_datasets(paths)
  df <- apply_cutoff(df, RSC)
  validate_df(df)

  if (flattened == "flattened") {
    df["NGS_read_count"] <- 1
  }

  df$del_length <- (df$End-1) - df$Start
  df$Shift <- df$del_length %% 3

  plot_df <- df %>%
    group_by(name, Shift) %>%
    summarise(counts=sum(NGS_read_count)) %>%
    mutate(Freq=counts / sum(counts) * 100) %>%
    mutate(Shift=case_when(
      Shift == 0 ~ "in-frame",
      Shift == 1 ~ "shift +1",
      Shift == 2 ~ "shift -1",
      TRUE ~ as.character(Shift)
    ))

  expected <- c(1/3, 1/3, 1/3)
  labels <- c()
  for (name in sort(unique(plot_df$name))) {
    observed_values <- plot_df$Freq[plot_df[["name"]] == name]
    r <- chisq.test(observed_values, p=expected)
    label <- paste(name, get_stat_symbol(r$p.value))
    labels <- c(labels, label)
  }

  pl <- ggplot(plot_df, aes(x=name, y=Freq, fill=factor(Shift))) +
    geom_bar(stat="identity", position="stack") +
    labs(x="Dataset" , y="Deletion shift [%]", fill="Shifts") +
    scale_x_discrete(labels=labels)

  prg$set(0.2, "Frame shifts plot")
  ggplotly(pl)
}


plot_multiple_segment_distribution <- function(paths, flattened, RSC, prg) {
  validate_selection(paths)
  df <- load_all_datasets(paths)
  df <- apply_cutoff(df, RSC)
  validate_df(df)

  if (flattened == "flattened") {
    df["NGS_read_count"] <- 1
  }

  plot_df <- df %>%
    group_by(name, strain, Segment) %>%
    summarise(counts=sum(NGS_read_count)) %>%
    mutate(Freq=counts / sum(counts) * 100) %>%
    rowwise() %>%
    mutate(seq_len=get_seq_len(strain, Segment))

  labels <- c()
  for (name in sort(unique(plot_df$name))) {
    subset_df <- plot_df[plot_df[["name"]] == name, ]
    expected <- subset_df$seq_len / sum(subset_df$seq_len)
    r <- chisq.test(subset_df$Freq, p=expected)
    label <- paste(name, "", get_stat_symbol(r$p.value))
    labels <- c(labels, label)
  }

  pl <- ggplot(plot_df, aes(x=name, y=Freq, fill=Segment)) +
    geom_bar(stat="identity", position="stack") +
    labs(x="Dataset" , y="Segment of DelVG [%]", fill="Segment") +
    scale_x_discrete(labels=labels)

  prg$set(0.3, "Segment distribution plot")
  ggplotly(pl)
}


plot_multiple_deletion_length<-function(paths,segment,flattened,n_bins, RSC, prg) {
  validate_selection(paths)
  df <- load_all_datasets(paths)
  df <- apply_cutoff(df, RSC)
  validate_df(df)
  
  # slice df by segment, reformat and bind on position and NGS count
  df <- df[df$Segment == segment, ]
  df$seq_len <- apply(df,
    1,
    function(row) get_seq_len(row["strain"], row["Segment"])
  )
  df["Length"] <- df["Start"] + (df["seq_len"] - df["End"] + 1)

  # multiply each column by NGS count if data is unflattened
  if (flattened != "flattened") {
    df <- data.frame(lapply(df, rep, df$NGS_read_count))
  }

  validate_plotting(df, segment)
  
  pl <- ggplot(df, aes(x=Length, fill=name)) +
    geom_histogram(position="identity", alpha=0.3, bins=n_bins) +
    labs(x="Length of DI candidate", y="Number of occurrences", fill="Dataset")
  
  prg$set(0.5, "Deletion length plot")
  ggplotly(pl)
}


plot_multiple_nucleotide_enrichment<-function(paths,segment,pos,flat,nuc,RSC, prg) {
  validate_selection(paths)
  df <- load_all_datasets(paths)
  df <- apply_cutoff(df, RSC)
  df <- df[df$Segment == segment, ]
  validate_df(df)

  exp_df <- load_expected_data(paths)
  exp_df <- exp_df[exp_df$Segment == segment, ]

  position <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  diff <- c()
  p_vals <- c()
  symb_ys <- c()
  symb_y <- 1
  unique_names <- sort(unique(df$name))
  for (name in unique_names) {
    n_df <- df[df$name == name, ]
    exp_n_df <- exp_df[exp_df$name == name, ]

    if (nrow(n_df) == 0 || nrow(exp_n_df) == 0) {
      diff <- c(diff, rep(0, 10))
      next
    }
    
    strain <- sort(unique(n_df$strain))
    counts <- prepare_nuc_enr_data(n_df, segment, flat, strain, pos, nuc)
    exp_counts<-prepare_nuc_enr_data(exp_n_df, segment, flat, strain, pos, nuc)

    comb <- cbind(counts, exp_counts$rel_occurrence)
    current_names <- colnames(comb)
    current_names[length(current_names)] <- "rel_occ2"
    colnames(comb) <- current_names
    
    comb$diff <- comb$rel_occurrence - comb$rel_occ2
    diff <- c(diff, comb$diff)

    # calculate ANOVA for each position
    n <- min(nrow(n_df), 1000)
    for (i in position) {
      obs_nucs <- as.integer(counts[i, "rel_occurrence"] * n)
      exp_nucs <- as.integer(exp_counts[i, "rel_occurrence"] * n)
      nucs <- c(rep(1, times=obs_nucs), rep(0, times=n-obs_nucs),
               rep(1, times=exp_nucs), rep(0, times=n-exp_nucs))
      group <- c(rep("obs", n), rep("exp", n))
      dat <- data.frame(nucs=nucs, group=group)
      anova <- aov(nucs ~ group, data=dat)
      p_vals <- c(p_vals, summary(anova)[[1]][["Pr(>F)"]][1])

    }
    symb_ys <- c(symb_ys, rep(symb_y, each=10))
    symb_y <- symb_y + 1
  }

  symbols <- gsub("ns.", "", lapply(p_vals, get_stat_symbol))
  positions <- c(rep(1:10, length(unique_names)))
  dataset <- rep(unique_names, each=10)
  plot_df <- data.frame(position=positions, dataset=dataset, diff=diff)
   
  y_max <- length(unique_names) + 0.5
  x1 <- ifelse(pos == "Start", 8, 3)
  x2 <- ifelse(pos == "Start", 3, 8)
  x_min <- ifelse(pos == "Start", 5.5, 0)
  x_max <- ifelse(pos == "Start", 11, 5.5)
  labels <- c("5", "4", "3", "2", "1", "-1", "-2", "-3", "-4", "-5")
  if (pos == "End") {
    labels <- rev(labels)
  }

  pl <- ggplot(plot_df, aes(x=position, y=dataset, fill=diff)) +
    geom_tile() +
    scale_fill_gradient2(low="blue", mid="white", high="red") +
    labs(x="Position", y="Dataset", fill="\u0394 (obs. - exp.)") +
    scale_x_continuous(breaks=position, labels=labels) +
    scale_y_discrete(expand = expansion(add = c(0.5, 1)))  +
    annotate("text", x=x1, y=y_max+0.05, label="deleted sequence") +
    annotate("text", x=x2, y=y_max+0.05, label="remaining sequence") +
    annotate("text", x=positions, y=symb_ys, label=symbols)

  value <- ifelse(pos == "Start", 0.7, 0.9)
  prg$set(value, "Nucleotide enrichment plot")
  ggplotly(pl)
}


plot_multiple_direct_repeat<-function(paths, segment, flattened, RSC, prg) {
  validate_selection(paths)
  df <- load_all_datasets(paths)
  df <- apply_cutoff(df, RSC)
  df <- df[df$Segment == segment, ]
  validate_df(df)

  exp_df <- load_expected_data(paths)
  exp_df <- exp_df[exp_df$Segment == segment, ]

  unique_names <- sort(unique(df$name))
  position <- c(rep(0:6, length(unique_names)))
  dataset <- c()
  diff <- c()
  for (name in unique_names) {
    r_df <- df[df$name == name, ]

    if (flattened == "flattened") {
      r_df$NGS_read_count[r_df$NGS_read_count != 1] <- 1
      n_df <- r_df
    } else {
      n_df<-data.frame(Start=integer(),End=integer(),NGS_read_count=integer())
      for (i in 1:nrow(r_df)) {
        n_df <- rbind(n_df, r_df[rep(i, r_df[i,3]),])
      }
    }
    
    exp_n_df <- exp_df[exp_df$name == name, ]
    if (nrow(n_df) == 0 || nrow(exp_n_df) == 0) {
      diff <- c(diff, rep(0, 7))
      next
    }
    validate_plotting(n_df, segment)

    strain <- sort(unique(n_df$strain))
    seq <- get_seq(strain, segment)
    n_samples <- nrow(n_df)
    n_df$direct_repeats <- apply(n_df,
      1,
      direct_repeats_counting_routine,
      seq
    )
    exp_n_df$direct_repeats <- apply(exp_n_df,
      1,
      direct_repeats_counting_routine,
      seq
    )
    df_1 <- prepare_direct_repeat_plot_data(n_df, "observed")
    df_2 <- prepare_direct_repeat_plot_data(exp_n_df, "expected")

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
    df_1$freq <- as.numeric(df_1$freq)
    df_2$freq <- as.numeric(df_2$freq)
    diff <- c(diff, df_1$freq - df_2$freq)

    # statistical testing with Chi-squared test
    matrix <- matrix(c(df_1$freq, df_2$freq), ncol=2) * 100
    r <- chisq.test(matrix)
    s <- get_stat_symbol(r$p.value)
    dataset <- c(dataset, rep(paste (name, s), each=7))
  }

  plot_df <- data.frame(position=position, dataset=dataset, diff=diff)
  position <- c(0:6)
  labels <- c("0", "1", "2", "3", "4", "5", ">5")

  pl <- ggplot(plot_df, aes(x=position, y=dataset, fill=diff)) +
    geom_tile() +
    scale_fill_gradient2(low="blue", mid="white", high="red") +
    labs(x="Direct repeat length", y="Dataset", fill="\u0394 (obs. - exp.)") +
    scale_x_continuous(breaks=position, labels=labels)

  prg$close()
  ggplotly(pl)
}

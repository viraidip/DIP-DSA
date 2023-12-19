load_single_dataset <- function(fname, sep=",") {
  path <- file.path(DATASETSPATH, fname)
  names <- c("Segment", "Start", "End", "NGS_read_count")
  cl <- c("character", "integer", "integer", "integer")
  if (file.exists(path)) {
    df <- read.csv(path, na.strings=c("NaN"), col.names=names, colClasses=cl, sep=sep)
  } else {
    df <- data.frame(
      "Segment"=character(),
      "Start"=integer(),
      "End"=integer(),
      "NGS_read_count"=integer()
    )
  }
  return(df)
}

load_all_datasets <- function(paths, sep=",") {
  df_list <- list()
  for (p in paths) {
    df <- load_single_dataset(p, sep=sep)
    df$name <- str_sub(str_extract(p, "/(.*?)\\."), 2, -2)
    df$strain <- str_extract(p, ".*(?=/)")
    df_list[[length(df_list) + 1]] <- df
  }
  final_df <- do.call(rbind, df_list)

  return(final_df)
}

load_expected_data <- function(paths){
  paths_tsv <- gsub("\\.csv$", ".tsv", paths)
  df <- load_all_datasets(paths_tsv, sep="\t")

  return(df)
}

plot_multiple_ngs_distribution <- function(paths, RCS) {
  df <- load_all_datasets(paths)
  df <- apply_cutoff(df, RCS)

  pl <- ggplot(df, aes(x=name, y=NGS_read_count, fill=name)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(title = "Distribution of NGS counts",
       x = "Dataset",
       y = "NGS count (log scale)")

  pl
}

plot_multiple_segment_distribution <- function(paths, flattened, RCS) {
  df <- load_all_datasets(paths)
  df <- apply_cutoff(df, RCS)

  percentage_df <- prop.table(table(df$name, df$Segment), margin = 1) * 100
  percentage_df <- as.data.frame(percentage_df)
  percentage_df$name <- rownames(percentage_df)

  # Create a barplot with percentage distribution
  pl <- ggplot(percentage_df, aes(x=Var1, y=Freq, fill=Var2)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x="dataset" , y="Segment of DVG [%]") +
    ggtitle("Segment where DVGs are originating from") + 
    theme(plot.title = element_text(size=20))

  ggplotly(pl)
}

plot_multiple_deletion_shift <- function(paths, flattened, RCS) {
  df <- load_all_datasets(paths)
  df <- apply_cutoff(df, RCS)

  df$del_length <- (df$End-1) - df$Start
  df$shift <- df$del_length %% 3

  percentage_df <- prop.table(table(df$name, df$shift), margin = 1) * 100
  percentage_df <- as.data.frame(percentage_df)
  percentage_df$name <- rownames(percentage_df)

  pl <- ggplot(percentage_df, aes(x=Var1, y=Freq, fill=Var2)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x="dataset" , y="Segment of DVG [%]") +
    ggtitle("Segment where DVGs are originating from") + 
    theme(plot.title = element_text(size=20))

  ggplotly(pl)
}


plot_multiple_deletion_length<-function(paths,segment,flattened,n_bins, RCS) {
  df <- load_all_datasets(paths)
  df <- apply_cutoff(df, RCS)
  
  # slice df by segment, reformat and bind on position and NGS count
  df <- df[df$Segment == segment, ]
  df$seq_len <- apply(df, 1, function(row) get_seq_len(row["strain"], row["Segment"]))
  df["Length"] <- df["Start"] + (df["seq_len"] - df["End"] + 1)

  # multiply each column by NGS count if data is unflattened
  if (flattened != "flattened") {
    df <- data.frame(lapply(df, rep, df$NGS_read_count))
  }

  validate_plotting(df, segment)
  
  pl <- ggplot(df, aes(x=Length, fill=name)) +
    geom_histogram(position="identity", alpha=0.3, bins=n_bins) +
    xlab("Length of DI candidate") +
    ylab("Number of occurrences") +
    ggtitle(paste("Histogram of DI RNA candidate lengths for segment",
      segment
      )
    ) + 
    theme(plot.title = element_text(size=20))

  ggplotly(pl)
}


plot_multiple_nucleotide_enrichment<-function(paths, segment, pos, flattened, nuc, RCS) {
  df <- load_all_datasets(paths)
  exp_df <- load_expected_data(paths)

  df <- apply_cutoff(df, RCS)
  df <- df[df$Segment == segment, ]  
  exp_df <- exp_df[exp_df$Segment == segment, ]

  unique_names <- unique(df$name)

  position <- c(rep(1:10, length(unique_names)))
  dataset <- rep(unique_names, each = 10)
  diff <- c()

  for (name in unique_names) {
    n_df <- df[df$name == name, ]
    exp_n_df <- exp_df[exp_df$name == name, ]

    if (nrow(n_df) == 0 || nrow(exp_n_df) == 0) {
      diff <- c(diff, rep(0, 10))
      next
    }
    
    strain <- unique(n_df$strain)
    counts <- prepare_nuc_dist_data(n_df, segment, flattened, strain, pos, nuc)
    exp_counts <- prepare_nuc_dist_data(exp_n_df, segment, flattened, strain, pos, nuc)

    comb <- cbind(counts, exp_counts$rel_occurrence)
    current_names <- colnames(comb)
    current_names[length(current_names)] <- "rel_occ2"
    colnames(comb) <- current_names
    
    comb$diff <- comb$rel_occurrence - comb$rel_occ2
    diff <- c(diff, comb$diff)

  }

  plot_df <- data.frame(position=position, dataset=dataset, diff=diff)

  position <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  labels <- c("5", "4", "3", "2", "1", "-1", "-2", "-3", "-4", "-5")
  if (pos == "End") {
    labels <- rev(labels)
  }
   
  y_max <- length(unique_names) + 0.5
  x1 <- ifelse(pos == "Start", 8, 3)
  x2 <- ifelse(pos == "Start", 3, 8)
  x_min <- ifelse(pos == "Start", 5.5, 0)
  x_max <- ifelse(pos == "Start", 11, 5.5)

  pl <- ggplot(plot_df, aes(x=position, y=dataset, fill=diff)) +
    geom_tile() +
    scale_fill_gradient2(low="blue", mid="white", high="red") +
    labs(x = "Position",
        y = "Dataset") +
    theme_minimal() +
    scale_x_continuous(breaks=position, labels=labels) +
   # annotate("text", x=position, y=y_text, label=symbols) +
    geom_rect(xmin=x_min, xmax=x_max, ymin=-1, ymax=y_max+0.5, alpha=0.5, fill="grey")+
    annotate("text", x=x1, y=y_max, label="deleted sequence") +
    annotate("text", x=x2, y=y_max, label="remaining sequence")

  ggplotly(pl)
}





plot_multiple_direct_repeat<-function(paths, segment, flattened, RCS) {
  df <- load_all_datasets(paths)
  exp_df <- load_expected_data(paths)

  df <- apply_cutoff(df, RCS)
  df <- df[df$Segment == segment, ]  
  exp_df <- exp_df[exp_df$Segment == segment, ]

  unique_names <- unique(df$name)

  position <- c(rep(0:6, length(unique_names)))
  dataset <- rep(unique_names, each=7)
  diff <- c()

  for (name in unique_names) {
    n_df <- df[df$name == name, ]
    exp_n_df <- exp_df[exp_df$name == name, ]

    if (nrow(n_df) == 0 || nrow(exp_n_df) == 0) {
      diff <- c(diff, rep(0, 7))
      next
    }
    
    strain <- unique(n_df$strain)
    seq <- get_seq(strain, segment)

    n_df["group"] <- rep("observed", nrow(n_df))
    exp_n_df["group"] <- rep("expected", nrow(exp_n_df))

    final_df <- rbind(n_df, exp_n_df)
    final_df["direct_repeats"] <- apply(final_df, 1, direct_repeats_counting_routine, seq)

    validate_plotting(final_df, segment)

    g1 <- unique(final_df[c("group")])[[1]][1]
    g2 <- unique(final_df[c("group")])[[1]][2]
    df_1 <- prepare_plot_data(final_df, g1)
    df_2 <- prepare_plot_data(final_df, g2)
    n_1 <- nrow(df_1)
    n_2 <- nrow(df_2)

    max_length <- max(max(df_1$length), max(df_2$length))
    for (i in 0:max_length) {
      if (!any(i==df_1[,1])) {
        df_1 <- rbind(df_1, list(i, 0.0, g1))
      }
      if (!any(i==df_2[,1])) {
        df_2 <- rbind(df_2, list(i, 0.0, g2))
      }
    }


    df_1$freq <- as.numeric(df_1$freq)
    df_2$freq <- as.numeric(df_2$freq)
    
    diff <- c(diff, df_1$freq - df_2$freq)

  }

  plot_df <- data.frame(position=position, dataset=dataset, diff=diff)
  position <- c(0:6)
  labels <- c("0", "1", "2", "3", "4", "5", ">5")

  pl <- ggplot(plot_df, aes(x=position, y=dataset, fill=diff)) +
    geom_tile() +
    scale_fill_gradient2(low="blue", mid="white", high="red") +
    labs(x="Direct repeat length",
        y="Dataset") +
    theme_minimal() +
    scale_x_continuous(breaks=position, labels=labels)

  ggplotly(pl)
}

get_intersecting_candidates <- function(df) {
  candidates <- c()
  for (name in unique(df$name)) {
    n_df <- df[df$name == name, ]
    cand <- paste(n_df$Segment, n_df$Start, n_df$End, sep="_")
    candidates <- append(candidates, cand)
  }

  final_df <- data.frame(table(candidates))
  colnames(final_df) <- c("DelVG candidate", "occurrence")
  thresh <- ceiling(length(unique(df$name)) / 2)
  final_df <- final_df[final_df$occurrence >= thresh, ]

  text <- paste("No candidates found in at least", thresh, "of the datasets.")
  shiny::validate(need((nrow(final_df) != 0), text)) 

  return(final_df)
}


plot_overlap_matrix <- function(paths, prg) {
  validate_selection(paths)
  df <- load_all_datasets(paths)
  validate_df(df)

  lists <- c()
  labels <- unique(df$name)
  for (name in labels) {
    n_df <- df[df$name == name, ]
    DI_list <- paste(n_df$Segment, n_df$Start, n_df$End, sep="_")
    lists <- c(lists, list(DI_list))
  }

  # initialize an empty matrix
  matrix_size <- length(lists)
  matrix <- matrix(0, nrow=matrix_size, ncol=matrix_size)
  # populate matrix entrywise
  for (i in 1:matrix_size) {
    set1 <- unique(lists[[i]])
    for (j in 1:matrix_size) {
      set2 <- unique(lists[[j]])
      matrix[i, j] <- (length(intersect(set1, set2)) / length(set2)) * 100
    }
  }

  plot_df <- as.data.frame(as.table(matrix))
  colnames(plot_df) <- c("Dataset1", "Dataset2", "intersection")
  p <- ggplot(plot_df, aes(x=Dataset1, y=Dataset2, fill=intersection)) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(x="Dataset", y="Dataset", fill="Intersecting DelVGs [%]") +
    scale_x_discrete(labels=labels) +
    scale_y_discrete(labels=labels)

  prg$set(0.3, "matrix plot")
  return(p)
}


plot_barplot_candidates <- function(paths, prg) {
  validate_selection(paths)
  df <- load_all_datasets(paths)
  validate_df(df)
  c_df <- get_intersecting_candidates(df)

  colnames(c_df) <- c("DelVG_candidate", "occurrence")
  segments <- str_split(c_df$DelVG_candidate, "_")
  c_df$segment <- sapply(segments, function(x) x[1])

  pl <- ggplot(c_df, aes(x=segment, fill=factor(occurrence))) +
    geom_bar(position="dodge") +
    labs(x="Segment", y="Identified DelVGs", fill="Datasets") +
    theme_minimal()
  
  prg$set(0.7, "candidate count plot")
  ggplotly(pl)
}


plot_highest_n_ranked <- function(paths, segment, thresh, prg) {
  validate_selection(paths)
  df <- load_all_datasets(paths)
  df <- df %>%
    group_by(name) %>%
    mutate(perc=ecdf(NGS_read_count)(NGS_read_count) * 100)

  c_df <- get_intersecting_candidates(df)
  df$key <- paste(df$Segment, df$Start, df$End, sep="_")
  marked_points <- df[df$key %in% as.character(c_df[["DelVG candidate"]]), ]
  marked_points <- marked_points[marked_points$Segment == segment, ]

  sum_ranking <- marked_points %>%
    group_by(key) %>%
    summarise(sum=sum(perc, na.rm=TRUE)) %>%
    arrange(desc(sum))

  mean_ranking <- marked_points %>%
    group_by(key) %>%
    summarise(mean=mean(perc, na.rm=TRUE)) %>%
    arrange(desc(mean))

  l1 <- sum_ranking$key
  l2 <- mean_ranking$key
  n_values <- c(1:min(thresh, length(l1)))
  overlap_counts <- numeric()
  found_cands <- c()
  cands_per_n <- c()
  for (n in n_values) {
    overlap <- intersect(l1[1:n], l2[1:n])
    overlap_counts <- c(overlap_counts, length(overlap))   
    new_cands <- setdiff(overlap, found_cands)
    if (length(new_cands) == 0) { # no new cand found
      new_cands <- " "
    } else if (length(new_cands) > 1) { # more than 1 new cand found
      found_cands <- c(found_cands, new_cands)
      new_cands <- paste(new_cands, collapse="\n")
    } else { # exactly 1 new cand found
      found_cands <- c(found_cands, new_cands)
    }    
    cands_per_n <- c(cands_per_n, new_cands)
  }

  results_df <- data.frame(n=n_values,
                           overlap_count=overlap_counts,
                           cands_per_n=cands_per_n)

  pl <- ggplot(results_df, aes(x=n, y=overlap_count)) +
    geom_line() +
    geom_point() +
    labs(x="n highest ranked DelVGs", y="# DelVGs in both rankings") +
    annotate("text",
             x=results_df$n+0.3,
             y=results_df$overlap_count-0.15,
             label=results_df$cands_per_n) +
    theme_minimal()

  prg$close()
  ggplotly(pl)
}

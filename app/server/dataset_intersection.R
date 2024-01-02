get_intersecting_candidates <- function(df) {
  candidates <- c()
  for (name in unique(df$name)) {
    n_df <- df[df$name == name, ]
    cand <- paste(n_df$Segment, n_df$Start, n_df$End, sep="_")
    candidates <- append(candidates, cand)
  }

  counts <- table(candidates)
  final_df <- data.frame(counts)
  colnames(final_df) <- c("DVG candidate", "occurrence")
  return(final_df)
}

intersecting_candidates_table <- function(paths, RSC, thresh) {
  df <- load_all_datasets(paths)
  df <- apply_cutoff(df, RSC)

  final_df <- get_intersecting_candidates(df)
  final_df <- final_df[final_df$occurrence >= thresh, ]
  return(final_df)
}


plot_overlap_matrix <- function(paths, RSC) {
  validation_text <- "No plot created. Please select at least two datasets."
  shiny::validate(need((length(paths) > 1), validation_text))

  df <- load_all_datasets(paths)
  df <- apply_cutoff(df, RSC)
  
  lists <- c()
  labels <- c()
  for (name in unique(df$name)) {
    n_df <- df[df$name == name, ]
    DI_list <- paste(n_df$Segment, n_df$Start, n_df$End, sep="_")
    lists <- c(lists, list(DI_list))
    labels <- c(labels, name)
  }

  # initialize an empty matrix
  matrix_size <- length(lists)
  matrix <- matrix(0, nrow = matrix_size, ncol = matrix_size)
  # populate matrix entrywise
  for (i in 1:matrix_size) {
    set1 <- unique(lists[[i]])
    for (j in 1:matrix_size) {
      set2 <- unique(lists[[j]])
      inter <- length(intersect(set1, set2))
      max <- max(length(set1), length(set2))
      matrix[i, j] <- (inter / max) * 100
    }
  }

  plot_df <- as.data.frame(as.table(matrix))
  colnames(plot_df) <- c("Dataset1", "Dataset2", "intersection")
  p <- ggplot(plot_df, aes(x=Dataset1, y=Dataset2, fill=intersection)) +
    geom_tile() +
    scale_fill_gradient(low="white", high="blue") +
    labs(x="Dataset", y="Dataset", fill="Intersection [%]") +
    scale_x_discrete(labels=labels) +
    scale_y_discrete(labels=labels)

  return(p)
}


plot_upset <- function(paths, RSC) {
  validation_text <- "No plot created. Please select at least two datasets."
  shiny::validate(need((length(paths) > 1), validation_text))

  df <- load_all_datasets(paths)
  df <- apply_cutoff(df, RSC)

  lists <- c()
  names <- c()
  for (name in unique(df$name)) {
    n_df <- df[df$name == name, ]
    DI_list <- paste(n_df$Segment, n_df$Start, n_df$End, sep="_")
    lists <- c(lists, list(DI_list))
    names <- c(names, name)
  }
  names(lists) <- names

  m = make_comb_mat(lists)
  p <- UpSet(m)
  return(p)
}


plot_intersecting_candidates_NGS <- function(paths, RSC) {
  df <- load_all_datasets(paths)
  df <- apply_cutoff(df, RSC)

  c_df <- get_intersecting_candidates(df)
  thresh <- max(c_df$occurrence)
  sel_df <- c_df[c_df$occurrence >= thresh, ]
  
  df$key <- paste(df$Segment, df$Start, df$End, sep="_")
  marked_points <- df[df$key %in% as.character(sel_df[["DVG candidate"]]), ]

  pl <- ggplot(df, aes(x=name, y=NGS_read_count, color=name)) +
    geom_boxplot() +
    scale_y_log10() +
    labs(x="Dataset", y="NGS count (log scale)", color="Legend", fill="") +
    geom_point(data=marked_points, aes(fill=key), size=3, shape=8)

  ggplotly(pl)
}

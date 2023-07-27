
load_dataset <- function(fname) {
  path <- file.path(DATASETSPATH, fname)
  names <- c("Segment", "Start", "End", "NGS_read_count")
  cl <- c("character", "integer", "integer", "integer")
  if (file.exists(path)) {
    df <- read.csv(path, na.strings=c("NaN"), col.names=names, colClasses=cl)
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

all_intersecting_candidates <- function(datasets, thresh) {
  candidates <- c()
  for (p in datasets) {
    df <- load_dataset(p)
    cand <- paste(df$Segment, df$Start, df$End, sep="_")
    candidates <- append(candidates, cand)
  }
  counts <- table(candidates)
  selected <- counts[counts >= thresh]
  return(data.frame(selected))
}

plot_overlap_matrix <- function(paths) {
  validation_text <- "No plot could be created. Please select at least two datasets."
  shiny::validate(need((length(paths) > 1), validation_text))
  
  lists <- c()
  for (p in paths) {
    df <- load_dataset(p)
    DI_list <- paste(df$Segment, df$Start, df$End, sep="_")
    lists <- c(lists, list(DI_list))
  }

  # initialize an empty matrix
  matrix_size <- length(lists)
  matrix <- matrix(0, nrow = matrix_size, ncol = matrix_size)

  # calculate the differences and populate the matrix
  for (i in 1:matrix_size) {
    for (j in 1:matrix_size) {
      set1 <- unique(lists[[i]])
      set2 <- unique(lists[[j]])

      matrix[i, j] <- length(intersect(set1, set2)) / (max(length(set1), length(set2)))
    }
  }
  
  split_paths <- strsplit(paths, "/")
  labels <- lapply(split_paths, function(path) {
    last_part <- tail(path, 1)
    cleaned_path <- gsub(".csv$", "", last_part)
    return(cleaned_path)
  })

  df <- as.data.frame(as.table(matrix))
  p <- ggplot(df, aes(x = Var1, y = Var2, fill=Freq)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(x = "Column", y = "Row", title = "Matrix Heatmap") +
    scale_x_discrete(labels=labels) +
    scale_y_discrete(labels=labels)

  return(p)
}

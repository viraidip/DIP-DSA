load_single_dataset <- function(fname) {
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

load_all_datasets <- function(paths) {
  df_list <- list()
  for (p in paths) {
    df <- load_single_dataset(p)
    name <- str_sub(str_extract(p, "/(.*?)\\."), 2, -2)
    df$name <- name
    df_list[[length(df_list) + 1]] <- df
  }
  
  final_df <- do.call(rbind, df_list)

  return(final_df)
}

load_expected_data <- function(paths){
  paths_tsv <- gsub("\\.csv$", ".tsv", paths)
  df <- load_all_datasets(paths_tsv)

  return(df)
}

plot_multiple_ngs_distribution <- function(paths, RCS) {
  df <- load_all_datasets(paths)
  df <- apply_cutoff(df, RCS)

  p <- ggplot(df, aes(x=name, y=NGS_read_count, fill=name)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(title = "Distribution of NGS counts",
       x = "Dataset",
       y = "NGS count (log scale)")

  p
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





### static values ###
SEGMENTS <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")

DATAPATH <- file.path(".", "data")
DATASETSPATH <- file.path(DATAPATH, "datasets")
FASTAPATH <- file.path(DATAPATH, "strain_segment_fastas")
TEMPPATH <- file.path(DATAPATH, "temp")

COLOR_MAP <- hash(A="blue", C="yellow", G="green", U="red")
NUC_MAP <- hash(A="Adenine", C="Cytosine", G="Guanine", U="Uracil")


### general functions ###
run_prechecks <- function() {
  # right now I have no prechecks implemented
}

get_seq <- function(strain, segment) {
  p <- file.path(DATASETSPATH,strain,"fastas", paste(segment,".fasta",sep=""))
  if (!file.exists(p)) {
    return(NULL)
  }
  fasta <- readDNAStringSet(p)
  return(RNAString(fasta[[1]]))
}

get_seq_len <- function(strain, segment) {
  return(length(get_seq(strain, segment)))
}

get_stat_symbol <- function(p) {
  if (is.nan(p)) {
    return("")
  } else if (p < 0.00001) {
    return("***")
  } else if (p < 0.001) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns.")
  }
}

packaging_signal_data_exists <- function(strain) {
  path <- file.path(DATASETSPATH, strain, "packaging_signal.json")
  return(file.exists(path))
}

load_packaging_signal_data <- function(strain) {
  path <- file.path(DATASETSPATH, strain, "packaging_signal.json")
  data <- read_json(path)
  return(data)
}

format_strain_name <- function(strain) {
# reformat name of strain folder
# needs to be called everytime the strain folder is accessed
# folder is e.g. A_California_07_2009
# but it is displayed and can be selected as A/California/07/2009
  return(gsub(pattern="/", replacement="_", x=strain))
}

validate_df <- function(df) {
  shiny::validate(need((nrow(df) != 0), "Empty dataframe")) 
}

validate_selection <- function(paths) {
  validation_text <- "No data selected. Please select at least one dataset."
  shiny::validate(need((length(paths) > 0), validation_text))
}

validate_plotting <- function(df, segment) {  
  validation_text <- paste(
    "No plot could be created. Segment",
    segment,
    "does not include any data points.",
    "Please select another one on the sidebar.")
  shiny::validate(need((nrow(df) != 0), validation_text))
}

apply_cutoff <- function(df, RCS) {
  df <- df[df$NGS_read_count >= RCS, ]
  return(df)
}

### random data generation ###
generate_sampling_data <- function(seq, s, e, n) {
  # create all combinations of start and end positions that are possible
  combs <- expand.grid(start=seq(s[1], s[2]), end=seq(e[1], e[2]))

  # create for each the DelVG Sequence
  sequences <- sapply(1:nrow(combs), function(i) {
    start_pos <- combs$start[i]
    end_pos <- combs$end[i]
    paste0(substr(seq, 1, start_pos), substr(seq, end_pos, nchar(seq)))
  })

  # create a data frame and select the indices where start is max for rows with
  # the same DelVG sequence
  temp_df <- data.frame(Start=combs$start, End=combs$end, Sequence=sequences)
  max_start_indices <- tapply(
    seq_len(nrow(temp_df)),
    temp_df$Sequence,
    function(indices) indices[which.max(temp_df$Start[indices])]
  )

  # create new dataframe with rows having max "Start" for each unique sequence
  # replace all other rows with the same sequence with this start/end
  # combination
  max_start_df <- temp_df[max_start_indices, ]
  max_start_df <- merge(temp_df, max_start_df, by="Sequence", all.x=TRUE)
  max_start_df <- max_start_df[, !duplicated(colnames(max_start_df))]

  # select only start and end column for sampling
  last_two_columns <- max_start_df[, -c(1:(ncol(max_start_df)-2))]
  names(last_two_columns) <- c("Start", "End")
  return(last_two_columns[sample(nrow(last_two_columns), n), , drop=FALSE])
}

create_random_data <- function(strain, dataset_name) {
  path <- file.path(DATASETSPATH, strain)
  file <- file.path(path, paste(dataset_name, ".csv", sep=""))
  names <- c("Segment", "Start", "End", "NGS_read_count")
  cl <- c("character", "integer", "integer", "integer")
  df <- read.csv(file, na.strings=c("NaN"), col.names=names, colClasses=cl)
  df <- apply_cutoff(df, 15) # filter here to remove at least some FP

  for (seg in SEGMENTS) {
    s_df <- df[df$Segment == seg, , drop=FALSE]
    if (nrow(s_df) == 0) {
      next
    }
    seq <- get_seq(strain, seg)
    start <- as.integer(mean(s_df$Start))
    end <- as.integer(mean(s_df$End))
    s <- c(max(start - 200, 50), start + 200)
    e <- c(end - 200, min(end + 200, nchar(seq) - 50))
    
    if (s[2] >= e[1] || s[1] == s[2] || e[1] == e[2]) {
      next
    }
    
    N_SAMPLES=5000
    if (exists("samp_df")) {
      temp_df <- generate_sampling_data(seq, s, e, N_SAMPLES)
      temp_df$Segment <- seg
      samp_df <- rbind(samp_df, temp_df)
    } else {
      samp_df <- generate_sampling_data(seq, s, e, N_SAMPLES)
      samp_df$Segment <- seg
    }
  }
  
  samp_df$NGS_read_count <- 1
  samp_df <- samp_df[, c(3, 1, 2, 4)]

  f_name <- paste(dataset_name, ".tsv", sep="")
  final_path <- file.path(path, f_name)
  write.table(samp_df, file=final_path, sep="\t", row.names=FALSE, quote=FALSE)
}


### data loading ###
load_single_dataset <- function(fname, sep=",") {
  path <- file.path(DATASETSPATH, fname)
  names <- c("Segment", "Start", "End", "NGS_read_count")
  cl <- c("character", "integer", "integer", "integer")
  if (file.exists(path)) {
    df <- read.csv(path,
      na.strings=c("NaN"),
      col.names=names,
      colClasses=cl,
      sep=sep
    )
  } else {
    df <- data.frame(
      "Segment"=character(),
      "Start"=integer(),
      "End"=integer(),
      "NGS_read_count"=integer()
    )
    validate_df(df)
  }
  return(df)
}

load_all_datasets <- function(paths, sep=",") {
  df_list <- list()
  for (p in paths) {
    df <- load_single_dataset(p, sep=sep)
    if (nrow(df) > 0){
      df$name <- str_sub(str_extract(p, "/(.*?)\\."), 2, -2)
      df$strain <- str_extract(p, ".*(?=/)")
      df_list[[length(df_list) + 1]] <- df      
    }
  }
  final_df <- do.call(rbind, df_list)

  return(final_df)
}

load_expected_data <- function(paths){
  paths_tsv <- gsub("\\.csv$", ".tsv", paths)
  df <- load_all_datasets(paths_tsv, sep="\t")

  return(df)
}


### single dataset tab ###
format_dataframe_lengths <- function(df, segment, strain, flattened) {
  df <- df[df$Segment == segment, ]
  seq_len <- get_seq_len(strain, segment)
  df["Length"] <- df["Start"] + (seq_len - df["End"] + 1)
  # multiply each column by NGS count if data is unflattened
  if (flattened != "flattened") {
    df <- data.frame(lapply(df, rep, df$NGS_read_count))
  }
  return(df)
}

add_stats_lengths <- function(df, pl) {
  col1 = "#8B2323"
  col2 = "#FF4040"
  y_f = 1.0
  
  # calculate stats and add them to plot
  mean <- mean(df$Length)
  median <- median(df$Length)
  mean_l <- paste("Mean =", format(mean, digits=5))
  median_l <- paste("Median =", format(median, digits=5))
  y <- max(ggplot_build(pl)$data[[1]]$count)
  pl <- pl +
    geom_vline(xintercept=mean, col=col1) +
    annotate("text", x=mean, y=y*y_f, label=mean_l, col=col1) +
    geom_vline(xintercept=median, col=col2) +
    annotate("text", x=median, y=y*(y_f-0.1), label=median_l, col=col2)
  return(pl)
}

format_dataframe_locations <- function(df, segment, flattened, strain) {
  # slice df by segment, reformat and bind on position and NGS count
  df <- df[df$Segment == segment, ]
  if (nrow(df) == 0) {
    return(data.frame())
  }

  start_df <- df[, c("Start", "NGS_read_count")]
  start_df["Class"] <- "Start"
  colnames(start_df)[colnames(start_df) == "Start"] <- "Position"
  end_df <- df[, c("End", "NGS_read_count")]
  end_df["Class"] <- "End"
  colnames(end_df)[colnames(end_df) == "End"] <- "Position"
  df <- rbind(start_df, end_df)

  # reformat if data should be displayed flattened
  if (flattened == "flattened") {
    df["NGS_read_count"] <- 1
  }

  df <- ddply(df,
    c("Position", "Class"),
    summarise,
    NGS_read_count=sum(NGS_read_count)
  )

  get_nuc_at_position <- function(position, seq) {
    as.character(subseq(seq, start=position, end=position))
  }
  sequence <- get_seq(strain, segment)
  df$Nucleotide <- sapply(df$Position, get_nuc_at_position, seq=sequence)

  return(df)
}

add_packaging_signal <- function(p, strain, segment) {
  packaging_signal <- load_packaging_signal_data(strain)
  slen <- get_seq_len(strain, segment)
  x <- unlist(packaging_signal[segment])
  y <- layer_scales(p)$y$get_limits()[2]
  color <- c("blue", "blue", "red", "red")
  p <- p + geom_vline(xintercept=x, color=color, linetype="dotted")
  return (p)
}


### direct repeats ###
direct_repeats_counting_routine <- function(row, sequence) {
  start <- as.integer(row["Start"])
  end <- as.integer(row["End"])
  i_s_1 <- start-5
  i_s_2 <- start
  i_e_1 <- end-5-1
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

prepare_direct_repeat_plot_data <- function(df, label) {
  table <- table(df$direct_repeats)
  df <- data.frame(table, rep(label, length(table)))
  colnames(df) <- c("length", "counts", "group")
  df$freq <- df$counts/sum(df$counts)
  df$length <- as.numeric(as.character(df$length))

  return(df)
}


### nucleotide enrichment ###
counting_routine <- function(l, window, letter, ngs_read_count) {
  count_indices <- unlist(gregexpr(letter, window))
  if (count_indices[1] != -1){
    indices <- unlist(count_indices)
    l[indices] <- l[indices] + ngs_read_count
  }
  return(l)
}

count_nuc_dist <- function(seq, positions, pos, ngs_read_counts, nuc) {
  count <- integer(10)
  for (i in 1:length(positions)) {
    p <- positions[[i]]
    ngs_read_count <- ngs_read_counts[[i]]
    if (pos == "Start") {
      window <- subseq(seq, start=p-4, end=p+5)  
    } else if (pos == "End") {
      window <- subseq(seq, start=p-5, end=p+4)
    }
    
    count <- counting_routine(count, window, nuc, ngs_read_count)
  }
  rel_occurrence <- count / sum(ngs_read_counts)
  position <- c(seq(1, 10))
  nucleotide <- c(rep(nuc, 10))
  return(data.frame(rel_occurrence, position, nucleotide))
}

prepare_nuc_enr_data <- function(df, segment, flattened, strain, pos, nuc) {
  df <- df[df$Segment == segment, ]
  if (nrow(df) == 0) {
    return(df)
  }
  
  if (flattened == "flattened") {
    df$NGS_read_count <- 1
  }
  ngs_read_counts <- df$NGS_read_count

  sequence <- get_seq(strain, segment)
  count_df <- count_nuc_dist(sequence, df[, pos], pos, ngs_read_counts, nuc)
  return(count_df)
}

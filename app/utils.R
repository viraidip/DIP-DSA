### DEFINING STATIC VALUES ###
SEGMENTS <- c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS")

DATAPATH <- file.path(".", "data")
DATASETSPATH <- file.path(DATAPATH, "datasets")
FASTAPATH <- file.path(DATAPATH, "strain_segment_fastas")
TEMPPATH <- file.path(DATAPATH, "temp")

#con <- file("../conda_path.txt","r")
#CONDAENV <- readLines(con,n=1)
#close(con)

COLOR_MAP <- hash(A="blue", C="yellow", G="green", U="red")
NUC_MAP <- hash(A="Adenine", C="Cytosine", G="Guanine", U="Uracil")


### DEFINING FUNCTIONS ###
run_prechecks <- function() {
  if (!file.exists(TEMPPATH)) {
    dir.create(TEMPPATH)
  }
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
  } else if (p < 0.0001) {
    return("****")
  } else if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("")
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

# reformat name of strain folder
# needs to be called everytime the strain folder is accessed
# folder is e.g. A_California_07_2009
# but it is displayed and can be selected as A/California/07/2009
format_strain_name <- function(strain) {
  return(gsub(pattern="/", replacement="_", x=strain))
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


generate_sampling_data <- function(seq, s, e, n) {
  # create all combinations of start and end positions that are possible
  combinations <- expand.grid(start=seq(s[1], s[2]), end=seq(e[1], e[2]))

  # create for each the DVG Sequence
  sequences <- sapply(1:nrow(combinations), function(i) {
    start_pos <- combinations$start[i]
    end_pos <- combinations$end[i]
    paste0(substr(seq, 1, start_pos-1), substr(seq, end_pos, nchar(seq)))
  })

  # create a data frame
  temp_df <- data.frame(Start=combinations$start, End=combinations$end, Sequence=sequences)
  max_start_indices <- tapply(
    seq_len(nrow(temp_df)),
    temp_df$Sequence,
    function(indices) indices[which.max(temp_df$Start[indices])]
  )

  # Create a new dataframe with rows having maximum "Start" for each unique "Sequence"
  max_start_df <- temp_df[max_start_indices, ]
  max_start_df <- merge(temp_df, max_start_df, by = "Sequence", all.x = TRUE)
  max_start_df <- max_start_df[, !duplicated(colnames(max_start_df))]
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

  for (seg in SEGMENTS) {
    temp_df <- df[df$Segment == seg, , drop = FALSE]
    if (nrow(temp_df) == 0) {
      next
    }
    seq <- get_seq(strain, seg)
    start <- as.integer(mean(temp_df$Start))
    end <- as.integer(mean(temp_df$End))
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
  
  f_name <- paste(dataset_name, ".tsv", sep="")
  final_path <- file.path(path, f_name)
  write.table(samp_df, file=final_path, sep="\t", row.names=FALSE)
}


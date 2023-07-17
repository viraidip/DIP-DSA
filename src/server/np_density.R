
create_np_plot <- function(df, strain, segment, areas) {
  # reformat data
  df <- format_dataframe_locations(df, segment, "flattened")
  a_df <- reformat_np_areas(areas)

  p <- ggplot(df, aes(x=Position, y=NGS_read_count, fill=Class)) +
    geom_bar(stat="identity", position="dodge", width=1) +
    xlim(0, get_seq_len(strain, segment)) +
    xlab("Nucleotide position on segment") +
    ylab("NGS read count") +
    ggtitle(paste("High NP density areas mapped to deletion sites for segment",
      segment
      )
    ) +
    theme(plot.title = element_text(size=20))

  # if available add high NP areas as rectangles to plot
  if (nrow(a_df) > 0) {
    p <- p + geom_rect(data=a_df,
      aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=label, alpha=0.9),
      inherit.aes=FALSE
    ) +
    scale_alpha(guide = "none")
  }

  ggplotly(p)
}

in_high_np_area <- function(row, a_df) {
  p <- as.integer(row["Position"])
  for (i in 1:nrow(a_df)) {
    s <- a_df[i, "start"]
    e <- a_df[i, "end"]
    if (p >= s && p <= e) {
      return("high")
    }
  }
  return("low")
}

count_area_occurrence <- function(df, a_df, label) {
  df["NP_area"] <- apply(df, 1, in_high_np_area, a_df=a_df)
  table_df <- data.frame(table(df$NP_area))
  # dividing low/all
  final_df <- data.frame(
    Label="Ratio",
    Freq=table_df[2, 2]/nrow(df),
    High=table_df[1, 2],
    Low=table_df[2, 2]
  )
  final_df["Class"] <- label
  return(final_df)
}

create_np_ratios_info <- function(df, strain, segment, areas) {
  # reformat data
  obs_df <- format_dataframe_locations(df, segment, "flattened")

  # create sampling data
  seq <- get_seq(strain, segment)
  n_samples <- nrow(df) * 5
  sam_df <- create_direct_repeat_sampling_data(df, n_samples, seq)
  sam_df["Segment"] <- rep(segment, n_samples)
  sam_df <- format_dataframe_locations(sam_df, segment, "flattened")

  # reformat np areas, if they are valid areas provided create a plot
  a_df <- reformat_np_areas(areas)
  
  if (nrow(a_df) == 0 || sum(is.na(a_df)) > 0) {
    return("")
  } else {
    # count occurrences inside and outside of NP areas
    o_df <- count_area_occurrence(obs_df, a_df, "observed")
    s_df <- count_area_occurrence(sam_df, a_df, "expected")
    r_df <- rbind(o_df, s_df)

    r1 <- r_df[1, "Freq"]
    r2 <- r_df[2, "Freq"]

    # statistical test (binom test)
    fill <- " not "
    if (!is.na(r_df[1, "Low"]) && !is.na(r2)) {
      pv <- binom.test(r_df[1, "Low"], nrow(obs_df), r2)$p.value
      pv <- round(pv, digits=8)
      if (pv < 0.05) {
        fill <- " "
      }
    } else {
      pv <- "NA"
    }

    return(
      paste(
        paste("ratio of observed data: ", round(r1, digits=2)),
        paste("ratio of expected data: ", round(r2, digits=2)),
        paste("p-value of binomial test: ", pv),
        paste("This means that there are", fill, "less deletion sites in low",
          " NP areas than expected.", sep=""),
        sep="\n"
      )
    )
  }

}


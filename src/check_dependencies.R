
pcks <- list("shiny", "shinydashboard", "ggplot2", "DT", "comprehenr", "hash", "plotly", "plyr", "jsonlite", "tools", "reticulate")

missing_pcks <- setdiff(pcks, rownames(installed.packages()))

for (m_pck in missing_pcks) {
  install.packages(m_pck)
}

if (!require("Biostrings")) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("Biostrings")
}


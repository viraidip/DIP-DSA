
pcks <- list("shiny",
  "shinydashboard",
  "ggplot2",
  "DT",
  "comprehenr",
  "hash",
  "plotly",
  "plyr",
  "jsonlite",
  "tools",
  "reticulate",
  "shinyvalidate",
  "ggvenn",
  "ComplexHeatmap",
  "reshape2"
)

missing_pcks <- setdiff(pcks, rownames(installed.packages()))

for (m_pck in missing_pcks) {
  install.packages(m_pck, repos="https://cloud.r-project.org")
}

if (!require("Biostrings")) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org")
  }
  BiocManager::install("Biostrings")
}


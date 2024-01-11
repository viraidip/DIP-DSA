pcks <- list("DT",
  "ggplot2",
  "hash",
  "jsonlite",
  "plotly",
  "plyr",
  "shiny",
  "shinydashboard",
  "shinyvalidate",
  "stringr",
  "dplyr"
)

missing_pcks <- setdiff(pcks, rownames(installed.packages()))
for (m_pck in missing_pcks) {
  install.packages(m_pck, repos="https://cloud.r-project.org")
}

# Bioconductor packages
if (!("BiocManager" %in% rownames(installed.packages()))) {
  install.packages("BiocManager", repos="https://cloud.r-project.org")
}
if (!("Biostrings" %in% rownames(installed.packages()))) {
  BiocManager::install("Biostrings")
}
if (!("ComplexHeatmap" %in% rownames(installed.packages()))) {
  BiocManager::install("ComplexHeatmap")
}

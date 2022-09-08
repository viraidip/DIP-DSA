
packages = list("shiny", "shinydashboard", "ggplot2", "DT", "shinyWidgets", "tools", "comprehenr",  "hash")

install.packages(setdiff(packages, rownames(installed.packages())))

if (!require("Biostrings")) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("Biostrings")
}

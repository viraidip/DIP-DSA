
run_prediction <- function(s, e, strain, segment, clf) {
  f_path <- file.path("py", "run_clf.py")
  source_python(f_path)
  label <- run_classification(s, e, strain, segment, clf)

  return(label)
}

show_clf_results <- function(label) {
  return(
    paste("Label:\t", label)
  )
}

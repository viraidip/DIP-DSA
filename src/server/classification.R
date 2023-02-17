
run_prediction <- function(s, e, strain, segment, clf) {
  f_path <- file.path("py", "run_clf.py")
  source_python(f_path)
  sequence <- toString(get_seq(format_strain_name(strain), segment))
  # reformat strain name
  label <- run_classification(s, e, strain, segment, sequence, clf)

  return(label)
}

show_clf_results <- function(label) {
  return(
    paste("Label:\t", label)
  )
}

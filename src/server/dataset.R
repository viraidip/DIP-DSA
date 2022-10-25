
move_files <- function(from_list, to_list) {
  for (i in 1:length(from_list)) {
    print(from_list[[i]])
    file.copy(
      from=from_list[[i]],
      to=to_list[[i]]
    )
  }
}

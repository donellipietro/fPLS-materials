## format

format_time <- function(t) {
  if (attr(t, "units") == "mins") {
    return(as.numeric(t) * 60)
  } else {
    return(as.numeric(t))
  }
}

## dataframes

add_results <- function(data, new, columns_names, groups_names = NULL) {
  n <- length(new[[1]])
  ## init
  if (is.null(groups_names)) {
    new_data <- data.frame(Group = paste(1:n))
  } else {
    new_data <- data.frame(Group = groups_names)
  }
  ## add columns
  for (subgroup in columns_names[-1]) {
    new_data <- cbind(new_data, new[[subgroup]])
  }
  colnames(new_data) <- columns_names
  ## append to given data
  data <- as.data.frame(rbind(data, new_data))
  colnames(data) <- columns_names
  return(data)
}

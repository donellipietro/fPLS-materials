## room for solutions
names_columns <- c("Group", names_approaches)
empty_df <- data.frame(matrix(NaN, nrow = 0, ncol = length(names_columns)))
colnames(empty_df) <- names_columns
times <- empty_df
rmses <- list()
rmses$centering <- empty_df
rmses$loadings <- empty_df
rmses$loadings_locs <- empty_df
rmses$scores <- empty_df
rmses$reconstruction <- empty_df
rmses$reconstruction_locs <- empty_df

## laod pers' results seqly
for (i in 1:n_reps) {
  ## laod batch and log if not present
  sink(file_log, append = TRUE)
  tryCatch(
    {
      path_batch <- paste(path_results, "batch_", i, "/", sep = "")
      load(paste(path_batch, "batch_", i, "_results_evaluation.RData", sep = ""))
    },
    error = function(e) {
      cat(paste("Error in test ", name_test, " - batch ", i, ": ", conditionMessage(e), "\n", sep = ""))
    }
  )
  sink()

  ## times
  times <- add_results(
    times,
    list(
      fPCA_off = format_time(results_evaluation$fPCA_off$execution_time),
      fPCA_gcv = format_time(results_evaluation$fPCA_gcv$execution_time),
      SpatialPCA = format_time(results_evaluation$SpatialPCA$execution_time)
    ),
    names_columns
  )
  ## rmse
  for (name in c("loadings", "loadings_locs", "scores", "reconstruction", "reconstruction_locs")) {
    rmses[[name]] <- add_results(
      rmses[[name]],
      list(
        fPCA_off = results_evaluation$fPCA_off$rmse[[name]],
        fPCA_gcv = results_evaluation$fPCA_gcv$rmse[[name]],
        SpatialPCA = results_evaluation$SpatialPCA$rmse[[name]]
      ),
      names_columns
    )
  }
  cat(paste("- Batch", i, "loaded\n"))
}

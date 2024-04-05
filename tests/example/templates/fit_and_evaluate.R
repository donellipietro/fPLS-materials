## room for solutions
results_evaluation <- list()

## load results_evaluation if available
if (file.exists(paste(path_batch, "batch_", i, "_results_evaluation.RData", sep = ""))) {
  load(paste(path_batch, "batch_", i, "_results_evaluation.RData", sep = ""))
}


### Model MV-PCA ----
model_MV_PCA <- NULL
file_model <- paste(path_batch, "batch_", i, "_fitted_model_MV_PCA.RData", sep = "")
if (file.exists(file_model) && !FORCE_FIT) {
  if(FORCE_EVALUATE) {
    cat("- Loading fitted MV-PCA ... \n")
    load(file_model)
  }
} else if ("MV_PCA" %in% names_models) {
  cat("- Fitting MV-PCA ... ")
  
  ## model fit
  start.time <- Sys.time()
  model_MV_PCA <- MV_PCA_wrapped(data, n_comp = n_comp)
  end.time <- Sys.time()
  cat(paste("finished after", end.time - start.time, attr(end.time - start.time, "units"), "\n"))
  
  ## add flags
  model_MV_PCA$model_traits$is_functional <- FALSE
  model_MV_PCA$model_traits$has_interpolator <- FALSE
  
  ## add execution time to results
  model_MV_PCA$results$execution_time <- end.time - start.time
  
  ## save fitted model
  save(
    index_batch = i,
    model_MV_PCA,
    file = file_model
  )
}
if(!is.null(model_MV_PCA)) {
  ## model evaluation
  results_evaluation$MV_PCA <- evaluate_results(model_MV_PCA, generated_data)
  
  ## clean the workspace
  rm(model_MV_PCA)
}


### Model fPCA ----
model_fPCA <- NULL
file_model <- paste(path_batch, "batch_", i, "_fitted_model_fPCA.RData", sep = "")
if (file.exists(file_model) && !FORCE_FIT) {
  if(FORCE_EVALUATE) {
    cat("- Loading fitted fPCA ... \n")
    load(file_model)
  }
} else if ("fPCA" %in% names_models) {
  cat("- Fitting fPCA ... ")
  
  ## model fit
  start.time <- Sys.time()
  model_fPCA <- fdaPDE2::fPCA(
    data = data,
    center = centering(calibrator = gcv(seed = seed, lambda = lambda_grid)),
    solver = sequential()
  )
  model_fPCA$fit(model_fPCA$results$calibration$lambda_centering_opt, n_pc = n_comp)
  end.time <- Sys.time()
  cat(paste("finished after", end.time - start.time, attr(end.time - start.time, "units"), "\n"))
  
  ## add flag
  model_fPCA$model_traits$is_functional <- TRUE
  model_fPCA$model_traits$has_interpolator <- TRUE
  
  ## add execution time to results
  model_fPCA$results$execution_time <- end.time - start.time
  
  ## save fitted model
  save(
    index_batch = i,
    model_fPCA,
    file = file_model
  )
  
}
if(!is.null(model_fPCA)) {
  ## evaluating fields at locations
  model_fPCA$results$loadings_locs <- model_fPCA$evaluate(model_fPCA$results$loadings)
  model_fPCA$results$X_mean_locs <- model_fPCA$evaluate(model_fPCA$results$X_mean)
  model_fPCA$results$X_hat_locs <- t(model_fPCA$evaluate(t(model_fPCA$results$X_hat)))
  
  ## model evaluation
  results_evaluation$fPCA <- evaluate_results(model_fPCA, generated_data)
  
  ## clean the workspace
  rm(model_fPCA)
}

### save results evaluation ----
save(
  ## batch index
  index_batch = i,
  ## results
  results_evaluation,
  ## path
  file = paste(path_batch, "batch_", i, "_results_evaluation.RData", sep = "")
)

## clean the workspace
if("generated_data" %in% ls()) rm(generated_data)
rm(results_evaluation)

cat(paste("- Batch", i, "compleated.\n"))

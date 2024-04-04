## room for solutions
results_evaluation <- list()

## load results_evaluation if available
if (file.exists(paste(path_batch, "batch_", i, "_results_evaluation.RData", sep = ""))) {
  load(paste(path_batch, "batch_", i, "_results_evaluation.RData", sep = ""))
}


### Molde MV-PCA ----
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
  model_MV_PCA <- MV_PCA_wrapped(data, center = FALSE, n_comp = n_comp)
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
  model_MV_PCA$results$X_mean_locs <- NULL
  results_evaluation$MV_PCA <- evaluate_results(model_MV_PCA, generated_data)
  
  ## clean the workspace
  rm(model_MV_PCA)
}


### Model fPCA no calibration ----
model_fPCA_off <- NULL
file_model <- paste(path_batch, "batch_", i, "_fitted_model_fPCA_off.RData", sep = "")
if (file.exists(file_model) && !FORCE_FIT) {
  if(FORCE_EVALUATE) {
    cat("- Loading fitted fPCA sequential (no calibration) ... \n")
    load(file_model)
  }
} else if ("fPCA_off" %in% names_models) {
  cat("- Fitting fPCA sequential (no calibration) ... ")
  
  ## model fit
  start.time <- Sys.time()
  model_fPCA_off <- fdaPDE2::fPCA(
    data = data,
    center = FALSE,
    solver = sequential()
  )
  model_fPCA_off$fit(hyperparameters(1e-7), n_pc = n_comp)
  end.time <- Sys.time()
  cat(paste("finished after", end.time - start.time, attr(end.time - start.time, "units"), "\n"))
  
  ## add flag
  model_fPCA_off$model_traits$is_functional <- TRUE
  model_fPCA_off$model_traits$has_interpolator <- TRUE
  
  ## add execution time to results
  model_fPCA_off$results$execution_time <- end.time - start.time
  
  ## save fitted model
  save(
    index_batch = i,
    model_fPCA_off,
    file = file_model
  )
  
}
if(!is.null(model_fPCA_off)) {
  ## evaluating fields at locations
  model_fPCA_off$results$loadings_locs <- model_fPCA_off$evaluate(model_fPCA_off$results$loadings)
  model_fPCA_off$results$X_mean_locs <- NULL
  model_fPCA_off$results$X_hat_locs <- t(model_fPCA_off$evaluate(t(model_fPCA_off$results$X_hat)))
  
  ## model evaluation
  results_evaluation$fPCA_off <- evaluate_results(model_fPCA_off, generated_data)
  
  ## clean the workspace
  rm(model_fPCA_off)
}


### Model fPCA (gcv calibration) ----
model_fPCA_gcv <- NULL
file_model <- paste(path_batch, "batch_", i, "_fitted_model_fPCA_gcv.RData", sep = "")
if (file.exists(file_model) && !FORCE_FIT) {
  if(FORCE_EVALUATE) {
    cat("- Loading fitted fPCA sequential (gcv calibration) ... \n")
    load(file_model)
  }
} else if ("fPCA_gcv" %in% names_models) {
  cat("- Fitting fPCA sequential (gcv calibration) ... ")
  
  ## model fit
  start.time <- Sys.time()
  model_fPCA_gcv <- fdaPDE2::fPCA(
    data = data,
    center = FALSE,
    solver = sequential()
  )
  model_fPCA_gcv$fit(calibrator = gcv(lambda = lambda_grid, seed = seed), n_pc = n_comp)
  end.time <- Sys.time()
  cat(paste("finished after", end.time - start.time, attr(end.time - start.time, "units"), "\n"))
  
  ## add flag
  model_fPCA_gcv$model_traits$is_functional <- TRUE
  model_fPCA_gcv$model_traits$has_interpolator <- TRUE
  
  ## add execution time to results
  model_fPCA_gcv$results$execution_time <- end.time - start.time
  
  ## save fitted model
  save(
    index_batch = i,
    model_fPCA_gcv,
    file = file_model
  )
}
if(!is.null(model_fPCA_gcv)) {
  ## evaluating fields at locations
  model_fPCA_gcv$results$loadings_locs <- model_fPCA_gcv$evaluate(model_fPCA_gcv$results$loadings)
  model_fPCA_gcv$results$X_mean_locs <- NULL
  model_fPCA_gcv$results$X_hat_locs <- t(model_fPCA_gcv$evaluate(t(model_fPCA_gcv$results$X_hat)))
  
  ## model evaluation
  results_evaluation$fPCA_gcv <- evaluate_results(model_fPCA_gcv, generated_data)
  
  ## clean the workspace
  rm(model_fPCA_gcv)
}


### Model fPCA monolithic (no calibration) ----
model_fPCA_mon_hyb <- NULL
file_model <- paste(path_batch, "batch_", i, "_fitted_model_fPCA_mon_hyb.RData", sep = "")
if (file.exists(file_model) && !FORCE_FIT) {
  if(FORCE_EVALUATE) {
    cat("- Loading fitted fPCA monolithic (hybrid calibration) ... \n")
    load(file_model)
  }
} else if ("fPCA_mon_hyb" %in% names_models) {
  cat("- Fitting fPCA monolithic (hybrid calibration) ... ")
  
  ## model fit
  start.time <- Sys.time()
  ## centering to find a suboptimal lambda (even if the data are already centered)
  model_centering <- fdaPDE2::fCentering(
    data = data,
    calibrator = gcv(seed = seed, lambda = lambda_grid)
  )
  model_centering$fit()
  ## real fit
  model_fPCA_mon_hyb <- fdaPDE2::fPCA(
    data = data,
    center = FALSE,
    solver = monolithic()
  )
  model_fPCA_mon_hyb$fit(model_centering$results$calibration$lambda_opt, n_pc = n_comp)
  end.time <- Sys.time()
  cat(paste("finished after", end.time - start.time, attr(end.time - start.time, "units"), "\n"))
  
  ## add flag
  model_fPCA_mon_hyb$model_traits$is_functional <- TRUE
  model_fPCA_mon_hyb$model_traits$has_interpolator <- TRUE
  
  ## add execution time to results
  model_fPCA_mon_hyb$results$execution_time <- end.time - start.time
  
  ## save fitted model
  save(
    index_batch = i,
    model_fPCA_mon_hyb,
    file = file_model
  )
}
if(!is.null(model_fPCA_mon_hyb)) {
  ## evaluating fields at locations
  model_fPCA_mon_hyb$results$loadings_locs <- model_fPCA_mon_hyb$evaluate(model_fPCA_mon_hyb$results$loadings)
  model_fPCA_mon_hyb$results$X_mean_locs <- NULL
  model_fPCA_mon_hyb$results$X_hat_locs <- t(model_fPCA_mon_hyb$evaluate(t(model_fPCA_mon_hyb$results$X_hat)))
  
  ## model evaluation
  results_evaluation$fPCA_mon_hyb <- evaluate_results(model_fPCA_mon_hyb, generated_data)
  
  ## clean the workspace
  rm(model_fPCA_mon_hyb)
}


### Model SpatialPCA ----
model_SpatialPCA <- NULL
file_model <- paste(path_batch, "batch_", i, "_fitted_model_SpatialPCA.RData", sep = "")
if (file.exists(file_model) && !FORCE_FIT) {
  if(FORCE_EVALUATE) {
    cat("- Loading fitted SpatialPCA ... \n")
    load(file_model)
  }
} else if ("SpatialPCA" %in% names_models) {
  cat("- Fitting SpatialPCA ... ")
  
  ## model fit
  model_SpatialPCA <- SpatialPCA_wrapped(data, n_comp)
  cat(paste("finished after", model_SpatialPCA$execution_time, attr(model_SpatialPCA$execution_time, "units"), "\n"))
  
  ## add flag
  model_SpatialPCA$model_traits$is_functional <- FALSE
  model_SpatialPCA$model_traits$has_interpolator <- TRUE
  
  ## add execution time to results
  model_SpatialPCA$results$execution_time <- model_SpatialPCA$execution_time
  
  ## save fitted model
  save(
    index_batch = i,
    model_SpatialPCA,
    file = file_model
  )
}
if(!is.null(model_SpatialPCA)) {
  ## model evaluation
  results_evaluation$SpatialPCA <- evaluate_results(model_SpatialPCA, generated_data)
  
  ## clean the workspace
  rm(model_SpatialPCA)
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
rm(generated_data, results_evaluation)

cat(paste("- Batch", i, "compleated.\n"))

## room for solutions
results_evaluation <- list()

### Molde 1 ----
cat("- Fitting MV-fPCA ... ")

## model fit
start.time <- Sys.time()
model_MV_PCA <- MV_PCA_wrapped(data, n_comp = n_comp)
end.time <- Sys.time()
cat(paste("finished after", end.time - start.time, attr(end.time - start.time, "units"), "\n"))

## model evaluation
results_evaluation$model1 <- evaluate_results(model_MV_PCA, generated_data)
results_evaluation$model1$execution_time <- end.time - start.time


### Molde 2 ----
cat("- Fitting fPCA (no calibration) ... ")

## model fit
start.time <- Sys.time()
model_fPCA_off <- fdaPDE2::fPCA(
  data = data,
  center = centering(calibrator = gcv(seed = seed, lambda = lambda_grid)),
  solver = sequential()
)
model_fPCA_off$fit(hyperparameters(1e-7), n_pc = n_comp)
end.time <- Sys.time()
cat(paste("finished after", end.time - start.time, attr(end.time - start.time, "units"), "\n"))

## evaluating fields at locations
model_fPCA_off$results$loadings_locs <- model_fPCA_off$evaluate(model_fPCA_off$results$loadings)
model_fPCA_off$results$X_mean_hat_locs <- t(model_fPCA_off$evaluate(t(model_fPCA_off$results$X_mean_hat)))
model_fPCA_off$results$X_hat_locs <- t(model_fPCA_off$evaluate(t(model_fPCA_off$results$X_hat)))

## model evaluation
results_evaluation$model2 <- evaluate_results(model_fPCA_off, generated_data)
results_evaluation$model2$execution_time <- end.time - start.time


### Molde 3 ----
cat("- Fitting fPCA (gcv calibration) ... ")

## model fit
start.time <- Sys.time()
model_fPCA_gcv <- fdaPDE2::fPCA(
  data = data,
  center = centering(calibrator = gcv(seed = seed, lambda = lambda_grid)),
  solver = sequential()
)
model_fPCA_gcv$fit(calibrator = gcv(lambda = lambda_grid, seed = seed), n_pc = n_comp)
end.time <- Sys.time()
cat(paste("finished after", end.time - start.time, attr(end.time - start.time, "units"), "\n"))

## evaluating fields at locations
model_fPCA_gcv$results$loadings_locs <- model_fPCA_gcv$evaluate(model_fPCA_gcv$results$loadings)
model_fPCA_gcv$results$X_mean_hat_locs <- t(model_fPCA_gcv$evaluate(t(model_fPCA_gcv$results$X_mean_hat)))
model_fPCA_gcv$results$X_hat_locs <- t(model_fPCA_gcv$evaluate(t(model_fPCA_gcv$results$X_hat)))


## model evaluation
results_evaluation$model3 <- evaluate_results(model_fPCA_gcv, generated_data)
results_evaluation$model3$execution_time <- end.time - start.time

## create batch directory
path_batch <- paste(path_results, "batch_", i, "/", sep = "")
mkdir(path_batch)

## save results evaluation
save(
  ## batch index
  index_batch = i,
  ## results
  results_evaluation,
  ## path
  file = paste(path_batch, "batch_", i, "_results_evaluation.RData", sep = "")
)

## save fitted models
save(
  ## batch index
  index_batch = i,
  ## fits
  model_MV_PCA, model_fPCA_off, model_fPCA_gcv,
  ## path
  file = paste(path_batch, "batch_", i, "_fitted_models.RData", sep = "")
)

## clean the worksapce
rm(generated_data, model_MV_PCA, model_fPCA_off, model_fPCA_gcv, results_evaluation)

cat(paste("- Batch", i, "compleated.\n"))

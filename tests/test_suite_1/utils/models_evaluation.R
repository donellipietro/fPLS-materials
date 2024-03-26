# performances evaluation ----


## utils

adjust_norms <- function(FF, SS) {
  n_comp <- ncol(SS)
  f_norms <- c()
  for (h in 1:n_comp) {
    f_norms[h] <- norm_l2(FF[, h])
    FF[, h] <- FF[, h] / f_norms[h]
    SS[, h] <- SS[, h] * f_norms[h]
  }
  return(list(loadings_locs = FF, scores = SS, norms = f_norms))
}


adjust_results <- function(F_hat_locs, S_hat, F_true_locs, Fs_hat_evaluated = NULL) {
  ## number of components
  n_comp <- ncol(S_hat)

  ## change signs to match the true ones
  for (h in 1:n_comp) {
    if (RMSE(F_hat_locs[, h] + F_true_locs[, h]) < RMSE(F_hat_locs[, h] - F_true_locs[, h])) {
      F_hat_locs[, h] <- -F_hat_locs[, h]
      S_hat[, h] <- -S_hat[, h]
      if (!is.null(Fs_hat_evaluated)) {
        for (i in 1:length(Fs_hat_evaluated)) {
          Fs_hat_evaluated[[i]][, h] <- -Fs_hat_evaluated[[i]][, h]
        }
      }
    }
  }

  ## normalize results at locations
  adjusted_results <- adjust_norms(F_hat_locs, S_hat)

  ## adjust data at nodes accordingly
  for (h in 1:n_comp) {
    if (!is.null(Fs_hat_evaluated)) {
      for (i in 1:length(Fs_hat_evaluated)) {
        Fs_hat_evaluated[[i]][, h] <- Fs_hat_evaluated[[i]][, h] / adjusted_results$norms[h]
      }
    }
  }

  return(list(
    loadings_locs = adjusted_results$loadings_locs,
    scores = adjusted_results$scores,
    loadings_evaluated_list = Fs_hat_evaluated
  ))
}

evaluate_results <- function(model, generated_data) {
  ## number of computed components
  n_comp <- generated_data$dimensions$n_comp

  ## room for results
  rmse <- list()
  irmse <- list()
  angles <- list()

  ## centering
  norm <- ifelse(RMSE(generated_data$X_mean_true_locs) == 0, 1, RMSE(generated_data$X_mean_true_locs))
  rmse$centering <- RMSE(model$results$X_mean_locs - generated_data$X_mean_true_locs) / norm
  # norm <- ifelse(RMSE(generated_data$X_mean_true_locs) == 0, 1, IRMSE(generated_data$X_mean_true, model))
  # irmse$centering <- ifelse(!is.null(model$results$X_mean) && !is.null(model$Psi()), IRMSE(model$results$X_mean - generated_data$X_mean_true, model) / norm, NULL)

  ## loadings & scores
  adjusted_results <- adjust_results(
    model$results$loadings_locs,
    model$results$scores,
    generated_data$loadings_true_locs,
    ifelse(!is.null(model$results$loadings), list(loadings = model$results$loadings), NULL)
  )
  for (h in 1:n_comp) {
    rmse$loadings_locs[h] <- RMSE(adjusted_results$loadings_locs[, h] - generated_data$loadings_true_locs[, h]) / RMSE(generated_data$loadings_true_locs[, h])
    rmse$scores[h] <- RMSE(adjusted_results$scores[, h] - generated_data$scores_true[, h]) / RMSE(generated_data$scores_true[, h])
    rmse$loadings[h] <- ifelse(!is.null(adjusted_results$loadings_evaluated_list$loadings), RMSE(adjusted_results$loadings_evaluated_list$loadings[, h] - generated_data$loadings_true[, h]), NULL)
    # irmse$loadings[h] <- ifelse(!is.null(adjusted_results$loadings_evaluated_list$loadings) && !is.null(model$Psi()), IRMSE(adjusted_results$loadings_evaluated_list$loadings[, h] - generated_data$loadings_true[, h], model), NULL)
  }

  ## data reconstruction
  rmse$reconstruction_locs <- RMSE(model$results$X_hat_locs - generated_data$X_true_locs) / RMSE(generated_data$X_true_locs)
  rmse$reconstruction <- ifelse(!is.null(model$results$X_hat), RMSE(model$results$X_hat - generated_data$X_true) / RMSE(generated_data$X_true), NULL)
  # irmse$reconstruction <- ifelse(!is.null(model$results$X_hat) && !is.null(model$Psi()), IRMSE(model$results$X_hat - generated_data$X_true, model) / IRMSE(generated_data$X_true, model), NULL)

  ## angles
  ## ....

  return(list(
    rmse = rmse,
    irmse = irmse,
    angles = angles
  ))
}

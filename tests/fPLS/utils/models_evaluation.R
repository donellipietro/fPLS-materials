# performances evaluation ----

## utils

adjust_norms <- function(X_space_directions_locs,
                         X_latent_scores,
                         X_loadings_locs, Y_loadings,
                         mode) {
  
  # X_latent_scores = X * Psi * X_space_directions
  # => X_latent_scores / norm = X * Psi * X_space_directions / norm
  # Y_latent_scores = Y * Y_space_directions
  # => this does not need re-normalization
  
  # X_hat = X_latent_scores * X_loadings^T
  # => X_hat = (X_latent_scores / norm) * (norm * X_loadings^T)
  
  # moreover, in PLS-R
  # Y_hat = X_latent_scores * Y_loadings^T
  # => Y_hat = (X_latent_scores / norm) * (norm * Y_loadings^T)
  
  n_comp <- ncol(X_space_directions_locs)
  norms_locs <- c()
  for (h in 1:n_comp) {
    norms_locs[h] <- norm_l2(X_space_directions_locs[, h])
    X_space_directions_locs[, h] <- X_space_directions_locs[, h] / norms_locs[h]
    X_latent_scores[, h] <- X_latent_scores[, h] / norms_locs[h]
    X_loadings_locs[, h] <- X_loadings_locs[, h] * norms_locs[h]
    if(mode == "PLS-R" || mode == "fPLS-R") {
      Y_loadings[, i] <- Y_loadings[, i] * norms_locs[h]
    }
  }
  
  return(list(
    X_space_directions_locs = X_space_directions_locs,
    X_latent_scores = X_latent_scores,
    X_loadings_locs = X_loadings_locs, 
    Y_loadings = Y_loadings,
    norms_locs = norms_locs
  ))
}


adjust_results <- function(Y_space_directions, Y_space_directions_true,
                           Y_latent_scores,
                           X_space_directions_locs, X_space_directions_true_locs,
                           X_latent_scores, X_loadings_locs, Y_loadings,
                           X_space_directions_evaluated_list, X_loadings_evaluated_list,
                           mode) {
  
  ## number of components
  n_comp <- ncol(X_latent_scores)
  
  ## normalize results at locations
  adjusted_results <- adjust_norms(X_space_directions_locs,
                                   X_latent_scores,
                                   X_loadings_locs, Y_loadings,
                                   mode)
  X_space_directions_locs <- adjusted_results$X_space_directions_locs
  X_latent_scores <- adjusted_results$X_latent_scores
  X_loadings_locs <- adjusted_results$X_loadings_locs
  Y_loadings <- adjusted_results$Y_loadings
  norms_locs <- adjusted_results$norms_locs
  
  ## adjust evaluated quantities accordingly
  for (h in 1:n_comp) {
    if (!is.null(X_space_directions_evaluated_list)) {
      for (i in 1:length(X_space_directions_evaluated_list)) {
        X_space_directions_evaluated_list[[i]][, h] <- X_space_directions_evaluated_list[[i]][, h] / norms_locs[h]
      }
    }
    if (!is.null(X_loadings_evaluated_list)) {
      for (i in 1:length(X_loadings_evaluated_list)) {
        X_loadings_evaluated_list[[i]][, h] <- X_loadings_evaluated_list[[i]][, h] * norms_locs[h]
      }
    }
  }
  
  ## change signs to match the true ones
  for (h in 1:n_comp) {
    ## Y space
    if (RMSE(Y_space_directions[, h] + Y_space_directions_true[, h]) < RMSE(Y_space_directions[, h] - Y_space_directions_true[, h])) {
      Y_space_directions[, h] <- - Y_space_directions[, h]
      Y_latent_scores[, h] <- - Y_latent_scores[, h]
      if(!(mode == "PLS-R" || mode == "fPLS-R")) {
        Y_loadings[, i] <- - Y_loadings[, i]
      }
    }
    ## X space
    if (RMSE(X_space_directions_locs[, h] + X_space_directions_true_locs[, h]) < RMSE(X_space_directions_locs[, h] - X_space_directions_true_locs[, h])) {
      X_space_directions_locs[, h] <- - X_space_directions_locs[, h]
      X_latent_scores[, h] <- - X_latent_scores[, h]
      X_loadings_locs[, h] <- - X_loadings_locs[, h]
      if(mode == "PLS-R" || mode == "fPLS-R") {
        Y_loadings[, i] <- - Y_loadings[, i]
      }
      if (!is.null(X_space_directions_evaluated_list)) {
        for (i in 1:length(X_space_directions_evaluated_list)) {
          X_space_directions_evaluated_list[[i]][, h] <- - X_space_directions_evaluated_list[[i]][, h]
        }
      }
      if (!is.null(X_loadings_evaluated_list)) {
        for (i in 1:length(X_loadings_evaluated_list)) {
          X_loadings_evaluated_list[[i]][, h] <- - X_loadings_evaluated_list[[i]][, h]
        }
      }
    }
  }
  
  return(list(
    Y_space_directions = Y_space_directions,
    Y_latent_scores = Y_latent_scores,
    X_space_directions_locs = X_space_directions_locs,
    X_latent_scores = X_latent_scores,
    X_loadings_locs = X_loadings_locs,
    Y_loadings = Y_loadings,
    X_space_directions_evaluated_list = X_space_directions_evaluated_list,
    X_loadings_evaluated_list = X_loadings_evaluated_list
  ))
}


evaluate_fPLS_at_locations <- function(model, n_comp, mean) {
  # mean
  if(mean) {
    model$results$X_mean_locs <- model$evaluate(model_fPLS_R_off$results$X_mean)
  } else {
    model$results$Y_mean <- NULL
    model$results$X_mean_locs <- NULL
  }
  # reconstruction
  model$results$X_hat_locs <- list()
  for(h in 1:n_comp) {
    model$results$X_hat_locs[[h]] <- t(model$evaluate(t(model$results$X_hat[[h]])))
  }
  # Beta
  model$results$Beta_hat_locs <- list()
  for(h in 1:n_comp) {
    model$results$Beta_hat_locs[[h]] <- model$evaluate(model$results$Beta_hat[[h]])
  }
  # directions & loadings
  model$results$X_space_directions_locs <- model$evaluate(model$results$X_space_directions)
  model$results$X_loadings_locs <- model$evaluate(model$results$X_loadings)
  
  return(model)
}


evaluate_results <- function(model, generated_data) {
  
  ## number of computed components
  n_comp <- generated_data$dimensions$n_comp
  
  ## room for results
  rmse <- list()
  irmse <- list()
  
  ## execution time
  execution_time <- model$results$execution_time
  
  ## RMSE Y and X (at locations) ----
  
  ## means
  if(!is.null(model$results$Y_mean)){
    norm <- ifelse(RMSE(generated_data$expected_results$Y_mean) == 0, 1, RMSE(generated_data$expected_results$Y_mean))
    rmse$Y_mean <- RMSE(model$results$Y_mean - generated_data$expected_results$Y_mean) / norm
  }
  if(!is.null(model$results$X_mean_locs)){
    norm <- ifelse(RMSE(generated_data$expected_results$X_mean_true_locs) == 0, 1, RMSE(generated_data$expected_results$X_mean_true_locs))
    rmse$X_mean_locs <- RMSE(model$results$X_mean_locs - generated_data$expected_results$X_mean_true_locs) / norm
  }
  
  ## data reconstruction
  for(h in 1:n_comp) {
    norm <- RMSE(generated_data$expected_results$Y_true)
    rmse$Y_reconstruction[h] <- RMSE(model$results$Y_hat[[h]] - generated_data$expected_results$Y_true) / norm
    norm <- RMSE(generated_data$expected_results$X_true_locs)
    rmse$X_reconstruction_locs[h] <- RMSE(model$results$X_hat_locs[[h]] - generated_data$expected_results$X_true_locs) / norm
  }
  
  ## Beta
  if(model$MODE == "fPLS-R" || model$MODE == "PLS-R") {
    for(h in 1:n_comp) {
      norm <- RMSE(generated_data$expected_results$Beta_locs)
      rmse$Beta_locs[h] <- RMSE(model$results$Beta_hat_locs[[h]] - generated_data$expected_results$Beta_locs) / norm
    }
  }
  
  # plot.field_tile(locations, model$results$Beta_hat_locs[[1]][, 1], LEGEND = TRUE)
  # plot.field_tile(locations, model$results$Beta_hat_locs[[2]][, 1], LEGEND = TRUE)
  # plot.field_tile(locations, model$results$Beta_hat_locs[[3]][, 1], LEGEND = TRUE)
  # plot.field_tile(locations, model$results$Beta_hat_locs[[4]][, 1], LEGEND = TRUE)
  # plot.field_tile(locations, generated_data$expected_results$Beta_locs[, 1], LEGEND = TRUE) + ggtitle("true")
  # 
  # head(model$results$Beta_hat_locs[[4]])
  # head(generated_data$expected_results$Beta_locs)
  
  ## directions, loadings & scores
  if(model$model_traits$has_interpolator) {
    X_space_directions_evaluated_list <- list(X_space_directions = model$results$X_space_directions)
    X_loadings_evaluated_list <- list(X_loadings = model$results$X_loadings)
  }else {
    X_space_directions_evaluated_list <- NULL
    X_loadings_evaluated_list <- NULL
  }
  adjusted_results <- adjust_results(
    model$results$Y_space_directions,
    generated_data$expected_results$Y_space_directions_true,
    model$results$Y_latent_scores,
    model$results$X_space_directions_locs,
    generated_data$expected_results$X_space_directions_true_locs,
    model$results$X_latent_scores,
    model$results$X_loadings_locs,
    model$results$Y_loadings,
    X_space_directions_evaluated_list,
    X_loadings_evaluated_list,
    model$MODE
  )
  for (h in 1:n_comp) {
    ## directions
    norm <- RMSE(generated_data$expected_results$Y_space_directions_true[, h])
    rmse$Y_space_directions[h] <- RMSE(adjusted_results$Y_space_directions[, h] - generated_data$expected_results$Y_space_directions_true[, h]) / norm
    norm <- RMSE(generated_data$expected_results$X_space_directions_true_locs[, h])
    rmse$X_space_directions_locs[h] <- RMSE(adjusted_results$X_space_directions_locs[, h] - generated_data$expected_results$X_space_directions_true_locs[, h]) / norm
    ## latent scores
    norm <- RMSE(generated_data$expected_results$Y_latent_scores_true[, h])
    rmse$Y_latent_scores[h] <- RMSE(adjusted_results$Y_latent_scores[, h] - generated_data$expected_results$Y_latent_scores_true[, h]) / norm
    norm <- RMSE(generated_data$expected_results$X_latent_scores_true[, h])
    rmse$X_latent_scores[h] <- RMSE(adjusted_results$X_latent_scores[, h] - generated_data$expected_results$X_latent_scores_true[, h]) / norm
    ## loadings
    norm <- RMSE(generated_data$expected_results$Y_loadings_true[, h])
    rmse$Y_loadings[h] <- RMSE(adjusted_results$Y_loadings[, h] - generated_data$expected_results$Y_loadings_true[, h]) / norm
    norm <- RMSE(generated_data$expected_results$X_loadings_true_locs[, h])
    rmse$X_loadings_locs[h] <- RMSE(adjusted_results$X_loadings_locs[, h] - generated_data$expected_results$X_loadings_true_locs[, h]) / norm
  }
  
  ## RMSE at nodes (if possible) ----
  
  if (model$model_traits$has_interpolator) {
    
    ## means
    if(!is.null(model$results$X_mean)){
      norm <- ifelse(RMSE(generated_data$expected_results$X_mean_true) == 0, 1, RMSE(generated_data$expected_results$X_mean_true))
      rmse$X_mean <- RMSE(model$results$X_mean - generated_data$expected_results$X_mean_true) / norm
    }
    
    ## data reconstruction
    for(h in 1:n_comp) {
      norm <- RMSE(generated_data$expected_results$X_true)
      rmse$X_reconstruction[h] <- RMSE(model$results$X_hat[[h]] - generated_data$expected_results$X_true) / norm
    }
    
    ## Beta
    if(model$MODE == "fPLS-R" || model$MODE == "PLS-R") {
      for(h in 1:n_comp) {
        norm <- RMSE(generated_data$expected_results$Beta)
        rmse$Beta[h] <- RMSE(model$results$Beta_hat[[h]] - generated_data$expected_results$Beta) / norm
      }
    }
    
    ## directions & loadings
    for (h in 1:n_comp) {
      ## directions
      norm <- RMSE(generated_data$expected_results$X_space_directions_true[, h])
      rmse$X_space_directions[h] <- RMSE(adjusted_results$X_space_directions_evaluated_list$X_space_directions[, h] - generated_data$expected_results$X_space_directions_true[, h]) / norm
      ## loadings
      norm <- RMSE(generated_data$expected_results$X_loadings_true[, h])
      rmse$X_loadings[h] <- RMSE(adjusted_results$X_loadings_evaluated_list$X_loadings[, h] - generated_data$expected_results$X_loadings_true[, h]) / norm
    }
    
  }
  
  # ## IRMSE (if possible)
  # if (model$model_traits$is_functional) {
  #   
  #   ## mean
  #   if(!is.null(model$results$X_mean)){
  #     norm <- IRMSE(generated_data$expected_results$X_mean_true, model)
  #     norm <- ifelse(norm == 0, 1, norm)
  #     irmse$centering <- IRMSE(model$results$X_mean - generated_data$expected_results$X_mean_true, model) / norm
  #   }
  #   
  #   ## reconstruction
  #   for (h in 1:n_comp) {
  #     norm <- IRMSE(t(generated_data$expected_results$X_true), model)
  #     irmse$X_reconstruction[h] <- IRMSE(t(model$results$X_hat[[h]] - generated_data$expected_results$X_true), model) / norm
  #   }
  #   
  #   ## Beta
  #   if(model$MODE == "fPLS-R" || model$MODE == "PLS-R") {
  #     for(h in 1:n_comp) {
  #       norm <- IRMSE(generated_data$expected_results$Beta, model)
  #       irmse$Beta[h] <- IRMSE(model$results$Beta_hat[[h]] - generated_data$expected_results$Beta, model) / norm
  #     }
  #   }
  #   
  #   ## directions & loadings
  #   for (h in 1:n_comp) {
  #     ## directions
  #     norm <- IRMSE(generated_data$expected_results$X_space_directions_true[, h], model)
  #     irmse$X_space_directions[h] <- IRMSE(adjusted_results$X_space_directions_evaluated_list$X_space_directions[, h] - generated_data$expected_results$X_space_directions_true[, h], model) / norm
  #     ## loadings
  #     norm <- IRMSE(generated_data$expected_results$X_loadings_true[, h], model)
  #     irmse$X_loadings[h] <- IRMSE(adjusted_results$X_loadings_evaluated_list$X_loadings[, h] - generated_data$expected_results$X_loadings_true[, h], model) / norm
  #   }
  #   
  # }
  
  return(list(
    execution_time = execution_time,
    rmse = rmse,
    irmse = irmse
  ))
}

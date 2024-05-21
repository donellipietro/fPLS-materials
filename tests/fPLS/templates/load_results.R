
cat("\nLoading results for quantitative analysis ...\n")

## room for solutions
names_columns <- c("Group", names_models)
empty_df <- data.frame(matrix(NaN, nrow = 0, ncol = length(names_columns)))
colnames(empty_df) <- names_columns
times <- empty_df
rmses <- list()
irmses <- list()
angles <- list()

## laod pers' results seqly
for (i in 1:n_reps) {
  ## laod batch and log if not present
  tryCatch(
    {
      path_batch <- paste(path_results, "batch_", i, "/", sep = "")
      load(paste(path_batch, "batch_", i, "_results_evaluation.RData", sep = ""))
    },
    error = function(e) {
      cat(paste("Error in test ", name_test, " - batch ", i, ": ", conditionMessage(e), "\n", sep = ""))
    }
  )
  
  ## times
  times <- add_results(
    times, extract_new_results(results_evaluation, names_models, "execution_time"),
    names_columns
  )
  ## rmse
  names <- c("Y_reconstruction", "X_reconstruction_locs", "Beta_locs",
             "Y_space_directions", "X_space_directions_locs",
             "Y_latent_scores", "X_latent_scores",
             "Y_loadings", "X_loadings_locs")
  for (name in names) {
    rmses[[name]] <- add_results(
      rmses[[name]], extract_new_results(results_evaluation, names_models, c("rmse", name)),
      names_columns
    )
  }
  # ## irmse
  # for (name in c("centering", "loadings", "reconstruction")) {
  #   irmses[[name]] <- add_results(
  #     irmses[[name]], extract_new_results(results_evaluation, names_models, c("irmse", name)),
  #     names_columns
  #   )
  # }
  # ## angles
  # for (name in c("subspaces_m", "components_m", "components_f", "orthogonality_m", "orthogonality_f")) {
  #   angles[[name]] <- add_results(
  #     angles[[name]], extract_new_results(results_evaluation, names_models, c("angles", name)),
  #     names_columns
  #   )
  # }
  cat(paste("- Batch", i, "loaded\n"))
}



if(RUN$qualitative_analysis) {
  cat("\nLoading results for qualitative analysis ...\n")
  
  ## room for solutions
  scores <- list()
  X_space_directions <- list()
  X_space_directions_locs <- list()
  X_space_directions_HR <- list()
  X_loadings <- list()
  X_loadings_locs <- list()
  X_loadings_HR <- list()
  
  ## load true data
  path_true_data <- paste("data/fPLS/", paste(mode_MV, name_mesh, "b", Beta_index, "comp", n_comp, sep = "_"), ".RData", sep = "")
  load(path_true_data)
  
  ## load batches
  for (i in 1:n_reps) {
    
    ## laod batch and log if not present
    tryCatch(
      {
        path_batch <- paste(path_results, "/", "batch_", i, "/", sep = "")
        for(name_model in names_models) {
          load(paste(path_batch, "batch_", i, "_fitted_model_", name_model, ".RData", sep = ""))
        }
      },
      error = function(e) {
        cat(paste("Error in test ", name_test, " - batch ", i, ": ", conditionMessage(e), "\n", sep = ""))
      }
    )


    ## adjust X_space_directions_true & X_loadings_true 
    X_space_directions_true_locs <- cbind(
      evaluate_field(locations, X_space_directions_true_grid[,1], mesh),
      evaluate_field(locations, X_space_directions_true_grid[,2], mesh),
      evaluate_field(locations, X_space_directions_true_grid[,3], mesh),
      evaluate_field(locations, X_space_directions_true_grid[,4], mesh)
    )
    X_space_directions_true <- cbind(
      evaluate_field(domain$nodes, X_space_directions_true_grid[,1], mesh),
      evaluate_field(domain$nodes, X_space_directions_true_grid[,2], mesh),
      evaluate_field(domain$nodes, X_space_directions_true_grid[,3], mesh),
      evaluate_field(domain$nodes, X_space_directions_true_grid[,4], mesh)
    )
    X_space_directions_true_HR <- cbind(
      evaluate_field(grid_HR, X_space_directions_true_grid[,1], mesh),
      evaluate_field(grid_HR, X_space_directions_true_grid[,2], mesh),
      evaluate_field(grid_HR, X_space_directions_true_grid[,3], mesh),
      evaluate_field(grid_HR, X_space_directions_true_grid[,4], mesh)
    )
    X_loadings_true_locs <- cbind(
      evaluate_field(locations, X_loadings_true_grid[,1], mesh),
      evaluate_field(locations, X_loadings_true_grid[,2], mesh),
      evaluate_field(locations, X_loadings_true_grid[,3], mesh),
      evaluate_field(locations, X_loadings_true_grid[,4], mesh)
    )
    X_loadings_true <- cbind(
      evaluate_field(domain$nodes, X_loadings_true_grid[,1], mesh),
      evaluate_field(domain$nodes, X_loadings_true_grid[,2], mesh),
      evaluate_field(domain$nodes, X_loadings_true_grid[,3], mesh),
      evaluate_field(domain$nodes, X_loadings_true_grid[,4], mesh)
    )
    X_loadings_true_HR <- cbind(
      evaluate_field(grid_HR, X_loadings_true_grid[,1], mesh),
      evaluate_field(grid_HR, X_loadings_true_grid[,2], mesh),
      evaluate_field(grid_HR, X_loadings_true_grid[,3], mesh),
      evaluate_field(grid_HR, X_loadings_true_grid[,4], mesh)
    )
    for(h in 1:ncol(X_space_directions_true_locs)) {
      norm_locs <- norm_l2(X_space_directions_true_locs[, h])
      ## directions
      X_space_directions_true_locs[, h] <- X_space_directions_true_locs[, h] / norm_locs
      X_space_directions_true[, h] <- X_space_directions_true[, h] / norm_locs
      X_space_directions_true_HR[, h] <- X_space_directions_true_HR[, h] / norm_locs
      ## loadings
      X_loadings_true_locs[, h] <- X_loadings_true_locs[, h] * norm_locs
      X_loadings_true[, h] <- X_loadings_true[, h] * norm_locs
      X_loadings_true_HR[, h] <- X_loadings_true_HR[, h] * norm_locs
    }
    
    ## directions, loadings & scores MV_PLS
    adjusted_results_MV_PLS <- adjust_results(
      model_MV_PLS$results$Y_space_directions,
      Y_space_directions_true,
      model_MV_PLS$results$Y_latent_scores,
      model_MV_PLS$results$X_space_directions_locs,
      X_space_directions_true_locs,
      model_MV_PLS$results$X_latent_scores,
      model_MV_PLS$results$X_loadings_locs,
      model_MV_PLS$results$Y_loadings,
      NULL, NULL,
      model_MV_PLS$MODE
    )
    
    ## ajust loadings model_fPLS_off
    model_fPLS_off <- evaluate_fPLS_at_locations(model_fPLS_gcv, n_comp, mean)
    model_fPLS_off_X_space_directions_HR <- cbind(
      evaluate_field(grid_HR, model_fPLS_off$results$X_space_directions[,1], mesh),
      evaluate_field(grid_HR, model_fPLS_off$results$X_space_directions[,2], mesh),
      evaluate_field(grid_HR, model_fPLS_off$results$X_space_directions[,3], mesh),
      evaluate_field(grid_HR, model_fPLS_off$results$X_space_directions[,4], mesh)
    )
    model_fPLS_off_X_loadings_HR <- cbind(
      evaluate_field(grid_HR, model_fPLS_off$results$X_loadings[,1], mesh),
      evaluate_field(grid_HR, model_fPLS_off$results$X_loadings[,2], mesh),
      evaluate_field(grid_HR, model_fPLS_off$results$X_loadings[,3], mesh),
      evaluate_field(grid_HR, model_fPLS_off$results$X_loadings[,4], mesh)
    )
    model_fPLS_off_X_space_directions_evaluated_list <- list(
      X_space_directions = model_fPLS_off$results$X_space_directions,
      X_space_directions_HR = model_fPLS_off_X_space_directions_HR
    )
    model_fPLS_off_X_loadings_evaluated_list <- list(
      X_loadings = model_fPLS_off$results$X_loadings,
      X_loadings_HR = model_fPLS_off_X_loadings_HR
    )
    adjusted_results_fPLS_off <- adjust_results(
      model_fPLS_off$results$Y_space_directions,
      Y_space_directions_true,
      model_fPLS_off$results$Y_latent_scores,
      model_fPLS_off$results$X_space_directions_locs,
      X_space_directions_true_locs,
      model_fPLS_off$results$X_latent_scores,
      model_fPLS_off$results$X_loadings_locs,
      model_fPLS_off$results$Y_loadings,
      model_fPLS_off_X_space_directions_evaluated_list,
      model_fPLS_off_X_loadings_evaluated_list,
      model_fPLS_off$MODE
    )
    
    ## ajust loadings model_fPLS_gcv
    model_fPLS_gcv <- evaluate_fPLS_at_locations(model_fPLS_gcv, n_comp, mean)
    model_fPLS_gcv_X_space_directions_HR <- cbind(
      evaluate_field(grid_HR, model_fPLS_gcv$results$X_space_directions[,1], mesh),
      evaluate_field(grid_HR, model_fPLS_gcv$results$X_space_directions[,2], mesh),
      evaluate_field(grid_HR, model_fPLS_gcv$results$X_space_directions[,3], mesh),
      evaluate_field(grid_HR, model_fPLS_gcv$results$X_space_directions[,4], mesh)
    )
    model_fPLS_gcv_X_loadings_HR <- cbind(
      evaluate_field(grid_HR, model_fPLS_gcv$results$X_loadings[,1], mesh),
      evaluate_field(grid_HR, model_fPLS_gcv$results$X_loadings[,2], mesh),
      evaluate_field(grid_HR, model_fPLS_gcv$results$X_loadings[,3], mesh),
      evaluate_field(grid_HR, model_fPLS_gcv$results$X_loadings[,4], mesh)
    )
    model_fPLS_gcv_X_space_directions_evaluated_list <- list(
      X_space_directions = model_fPLS_gcv$results$X_space_directions,
      X_space_directions_HR = model_fPLS_gcv_X_space_directions_HR
    )
    model_fPLS_gcv_X_loadings_evaluated_list <- list(
      X_loadings = model_fPLS_gcv$results$X_loadings,
      X_loadings_HR = model_fPLS_gcv_X_loadings_HR
    )
    adjusted_results_fPLS_gcv <- adjust_results(
      model_fPLS_gcv$results$Y_space_directions,
      Y_space_directions_true,
      model_fPLS_gcv$results$Y_latent_scores,
      model_fPLS_gcv$results$X_space_directions_locs,
      X_space_directions_true_locs,
      model_fPLS_gcv$results$X_latent_scores,
      model_fPLS_gcv$results$X_loadings_locs,
      model_fPLS_gcv$results$Y_loadings,
      model_fPLS_gcv_X_space_directions_evaluated_list,
      model_fPLS_gcv_X_loadings_evaluated_list,
      model_fPLS_gcv$MODE
    )
    
    ## save adjusted directions at locations
    X_space_directions_locs$MV_PLS[[i]] <- adjusted_results_MV_PLS$X_space_directions_locs
    X_space_directions_locs$fPLS_off[[i]] <- adjusted_results_fPLS_off$X_space_directions_locs
    X_space_directions_locs$fPLS_gcv[[i]] <- adjusted_results_fPLS_gcv$X_space_directions_locs
    ## save adjusted directions at nodes
    X_space_directions$fPLS_off[[i]] <- adjusted_results_fPLS_off$X_space_directions_evaluated_list$X_space_directions
    X_space_directions$fPLS_gcv[[i]] <- adjusted_results_fPLS_gcv$X_space_directions_evaluated_list$X_space_directions
    ## save adjusted directions at grid
    X_space_directions_HR$fPLS_off[[i]] <- adjusted_results_fPLS_off$X_space_directions_evaluated_list$X_space_directions_HR
    X_space_directions_HR$fPLS_gcv[[i]] <- adjusted_results_fPLS_gcv$X_space_directions_evaluated_list$X_space_directions_HR
    
    ## save adjusted loadings at locations
    X_loadings_locs$MV_PLS[[i]] <- adjusted_results_MV_PLS$X_loadings_locs
    X_loadings_locs$fPLS_off[[i]] <- adjusted_results_fPLS_off$loadings_locs
    X_loadings_locs$fPLS_gcv[[i]] <- adjusted_results_fPLS_gcv$loadings_locs
    ## save adjusted loadings at nodes
    X_loadings$fPLS_off[[i]] <- adjusted_results_fPLS_off$X_loadings_evaluated_list$X_loadings
    X_loadings$fPLS_gcv[[i]] <- adjusted_results_fPLS_gcv$X_loadings_evaluated_list$X_loadings
    ## save adjusted loadings at grid
    X_loadings_HR$fPLS_off[[i]] <- adjusted_results_fPLS_off$X_loadings_evaluated_list$X_loadings_HR
    X_loadings_HR$fPLS_gcv[[i]] <- adjusted_results_fPLS_gcv$X_loadings_evaluated_list$X_loadings_HR

    cat(paste("- Batch", i, "loaded\n"))
  }
  
}

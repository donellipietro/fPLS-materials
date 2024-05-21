
## data generator
generate_2D_data <- function(
    
  ## true data identificator
  mode,     
  name_mesh,
  Beta_index,      
  n_comp,
  
  ## simulation setting
  domain,
  locs = NULL,
  n_stat_units = 50,
  n_nodes = 30,
  NSR_X_last_comp = 0.5,
  NSR_Y = 0.5,
  seed = 0,
  VERBOSE = FALSE) {
  
  ## load true data
  path_true_data <- paste("data/fPLS/", paste(mode, name_mesh, "b", Beta_index, "comp", n_comp, sep = "_"), ".RData", sep = "")
  load(path_true_data)
  
  ## set defaults
  if (is.null(locs)) {
    locs <- domain$nodes
  }
  
  ## nodes
  nodes <- domain$nodes
  
  ## dimensions
  n_nodes <- nrow(nodes)
  n_locs <- nrow(locs)
  n_resp <- nrow(Y_loadings_true)
  
  # X --------------------------------------------------------------------------
  
  ## generating the loadings functions
  X_loadings_true <- matrix(0, nrow = n_nodes, ncol = n_comp)
  X_loadings_true_locs <- matrix(0, nrow = n_locs, ncol = n_comp)
  X_space_directions_true <- matrix(0, nrow = n_nodes, ncol = n_comp)
  X_space_directions_true_locs <- matrix(0, nrow = n_locs, ncol = n_comp)
  for (i in 1:n_comp) {
    ## space directions evaluation
    X_space_directions_true_locs[, i] <- evaluate_field(locs, X_space_directions_true_grid[, i], grid)
    norm_locs <- RMSE(X_space_directions_true_locs[, i])
    X_space_directions_true_locs[, i] <- X_space_directions_true_locs[, i] / norm_locs
    X_space_directions_true[, i] <- evaluate_field(nodes, X_space_directions_true_grid[, i], grid) / norm_locs
    
    ## loadings evaluation
    X_loadings_true_locs[, i] <- evaluate_field(locs, X_loadings_true_grid[, i], grid) * norm_locs
    X_loadings_true[, i] <- evaluate_field(nodes, X_loadings_true_grid[, i], grid) * norm_locs
    
    ## consequent normalizations
    sd_X_latent_scores[i] <- sd_X_latent_scores[i] / norm_locs
    if(mode == "PLS-R") {
      Y_loadings_true[, i] <- Y_loadings_true[, i] * norm_locs
    }
  }
  
  ## Beta evaluation
  if(mode == "PLS-R") {
    Beta_true <- matrix(0, nrow = n_nodes, ncol = n_resp)
    Beta_true_locs <- matrix(0, nrow = n_locs, ncol = n_resp)
    for(i in 1:n_resp) {
      Beta_true_locs[, i] <- evaluate_field(locs, Beta_true_grid[[n_comp]][, i], grid)
      Beta_true[, i] <- evaluate_field(nodes, Beta_true_grid[[n_comp]][, i], grid)
    }
  }
  
  # noise X --------------------------------------------------------------------
  
  ###### -------> NEW
  # Computing the latent components sd. as 
  ## sd_comp_i = sqrt(Var[score_i*loading_i]) = sqrt(Var[score_i] * semi_range(loading_i)^2) 
  ## sd_comp_i = sqrt(Var[score_i]) * semi_range(loading_i)
  ## sd_comp_i = sd_X_latent_scores_i * semi_range(loading_i)
  semi_range_X_loadings_true <- 0.5 * (apply(X_loadings_true, 2, max)-apply(X_loadings_true, 2, min)) 
  sd_comp <- sd_X_latent_scores * semi_range_X_loadings_true
  
  ## computing the noise sd. (sigma_noise) for X
  ## NSR_X_last_comp = Var[noise]/Var[last comp.]
  ## Var[last comp.] = Var[score_last * X_loadings_last] = sd_comp_last^2
  ## => Var[noise] = NSR_X_last_comp * Var[3rd comp.] = NSR_X_last_comp * sd_s_3^2
  sigma_noise_x <- sqrt(NSR_X_last_comp * min(sd_comp^2))
  NSR_X <- NSR_X_last_comp * min(sd_comp^2) / sum(sd_comp^2)
  
  if (VERBOSE) {
    cat("\n\n# Computing the noise sd. (sigma_noise) for X")
    cat("\nNSR_X_last_comp = Var[noise]/Var[3rd comp.]")
    cat("\nVar[3rd comp.] = Var[score_3 * X_loadings_3] = sd_s_3^2")
    cat("\n=> Var[noise] = NSR_X_last_comp * Var[3rd comp.] = NSR_X_last_comp * sd_s_3^2")
    cat("\n=> sigma_noise = sqrt(Var[noise]) = sqrt(NSR_X_last_comp * sd_s_3^2")
    cat("\n")
    cat(paste("\n- Desired Noise to Signal Ratio:", NSR_X_last_comp))
    cat(paste("\n- Effective noise sd:", sigma_noise_x))
    cat(paste("\n- Global NSR_X:", NSR_X))
    cat("\n")
  }
  
  # noise Y----------------------------------------------------------------------
  
  ## computing the noise sd. (sigma_noise) for Y
  sigma_noise_y <- sqrt(NSR_Y) * as.numeric(sqrt(Y_loadings_true^2 %*% sd_X_latent_scores^2))
  
  if (VERBOSE) {
    cat("\n\n# Computing the noise sd. (sigma_noise) for Y")
    cat(paste("\n- Desired Noise to Signal Ratio:", NSR_Y))
    cat(paste("\n- Effective noise sd:", sigma_noise_y))
    cat("\n")
  }
  
  
  # sampling -------------------------------------------------------------------
  
  # Sampling scores and noise from a multivariate random normal with 0 mean 
  # and suitable variances as computed before
  Sampled <- mvrnorm(n_stat_units, mu = rep(0,  n_comp + n_locs + n_resp),
                     diag(c(sd_X_latent_scores^2, rep(sigma_noise_x^2, n_locs), sigma_noise_y)))
  X_latent_scores_true <- scale(as.matrix(Sampled[, (1):(n_comp)], ncol = n_comp), scale = FALSE)
  EE <- scale(as.matrix(Sampled[, (n_comp + 1):(n_comp + n_locs)], ncol = n_locs), scale = FALSE)
  FF <- scale(as.matrix(Sampled[, (n_comp + n_locs + 1):(n_comp + n_locs + n_resp)], ncol = n_resp),  scale = FALSE)
  
  
  # data -----------------------------------------------------------------------
  
  ## generating X:
  X_c_true <- X_latent_scores_true %*% t(X_loadings_true)
  X_c_true_locs <- X_latent_scores_true %*% t(X_loadings_true_locs)
  
  ## generating the mean function
  X_mean_true <- 0 # mean_generator(nodes)
  X_mean_true_locs <- 0 # mean_generator(locs)
  # X_mean_true <- X_mean_true * semi_range_X_c_true
  # X_mean_true_locs <- X_mean_true_locs * semi_range_X_c_true
  
  ## adding the mean
  X_true <- X_c_true # + rep(1, n_stat_units) %*% t(X_mean_true)
  X_true_locs <- X_c_true_locs # + rep(1, n_stat_units) %*% t(X_mean_true_locs)
  
  ## generating Y:
  Y_c_true <- X_latent_scores_true %*% t(Y_loadings_true)
  Y_mean_true <- 0 # t(c(2, 5)[1:L])
  Y_true <- Y_c_true # + rep(1,N) %*% Y_mean_true
  
  ## adding the noise and the mean:
  X_locs <- X_true_locs + EE
  Y <- Y_true + FF
  
  ## data
  data <- list(
    X = X_locs,
    Y = Y
  )
  
  ## expected results:
  expected_results <- list(
    ## data
    X_true = X_true,
    X_true_locs = X_true_locs,
    Y_true = Y_true,
    X_c_true = X_c_true,
    X_c_true_locs = X_c_true_locs,
    Y_c_true = Y_c_true,
    X_mean_true = X_mean_true,
    X_mean_true_locs = X_mean_true_locs,
    Y_mean_true = Y_mean_true,
    ## space directions
    X_space_directions_true = X_space_directions_true,
    X_space_directions_true_locs = X_space_directions_true_locs,
    Y_space_directions_true = Y_space_directions_true,
    ## loadings
    X_loadings_true = X_loadings_true,
    X_loadings_true_locs = X_loadings_true_locs,
    Y_loadings_true = Y_loadings_true,
    ## latent scores
    X_latent_scores_true = X_latent_scores_true,
    Y_latent_scores_true = X_latent_scores_true
  )
  if(mode == "PLS-R") {
    expected_results$Beta <- Beta_true
    expected_results$Beta_locs <- Beta_true_locs
  }
  
  return(list(
    ## dimensions
    dimensions = list(
      n_stat_units = n_stat_units,
      n_comp = n_comp,
      n_nodes = n_nodes,
      n_locs = n_locs,
      n_resp = n_resp
    ),
    ## data
    data = data,
    ## expected results
    expected_results = expected_results,
    ## computed
    sd_X_latent_scores = sd_X_latent_scores,
    sigma_noise_x = sigma_noise_x,
    sigma_noise_y = sigma_noise_y,
    NSR_X, NSR_Y
  ))
}




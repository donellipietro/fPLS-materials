#Copiato da file pietro

PLS <- function(Y, X, n_comp = 4, center = TRUE, MODE = "PLS-R"){
  
  if (MODE != "PLS-A" & MODE != "PLS-SB" & MODE != "PLS-R") {
    stop("Error: MODE must be either 'PLS-A', 'PLS-SB', or 'PLS-R'")
  }
  
  #### Initialization  ####
  
  # Centering
  if(center) {
    X_c = scale(X, scale = FALSE)
    X_mean = attr(X_c, "scaled:center")
    Y_c <- scale(Y, scale = FALSE)
    Y_mean = attr(Y_c, "scaled:center")
  } else {
    X_c = X
    X_mean = rep(0, ncol(X))
    Y_c = Y
    Y_mean = rep(0, ncol(Y))
  }
  
  #Extracting dimensions
  N <- nrow(X_c) #N
  L <- ncol(Y_c) #L
  S <- ncol(X_c) #S
  
  #X and Y directions
  W <- matrix(0, nrow = S, ncol = n_comp)
  V <- matrix(0, nrow = L, ncol = n_comp)
  #X and Y scores/components
  T <- matrix(0, nrow = N, ncol = n_comp)
  U <- matrix(0, nrow = N, ncol = n_comp)
  #X and Y loadings
  P <- matrix(0, nrow = S, ncol = n_comp)
  Q <- matrix(0, nrow = L, ncol = n_comp)
  
  #B matrix
  B <- matrix(0, nrow = n_comp, ncol = n_comp)
  Beta <- list()
  
  # Mean values in matrix form
  X_MEAN = matrix(X_mean, nrow = N, ncol = S, byrow = TRUE)
  Y_MEAN = matrix(Y_mean, nrow = N, ncol = L,  byrow = TRUE)
  
  # Matrices for intermediate steps
  XX <- list()
  YY <- list()
  
  # X and Y residual matrices
  EE <- X_c
  FF <- Y_c
  XX[[1]] <- EE
  YY[[1]] <- FF
  
  ## room for results
  Y_hat <- list()
  X_hat <- list()
  Beta_hat <- list()
  
  
  #### Main loop  ####
  
  for(h in 1:n_comp){
    
    ### Step 1: covariance maximization     ### 
    
    #Weight extraction
    C.YX <- t(FF) %*% EE
    SVD <- svd(C.YX, nu = 1, nv = 1)
    W[,h] <- SVD$v  
    V[,h] <- SVD$u 
    
    #Score/component computation
    T[,h] <- EE %*% W[,h] 
    U[,h] <- FF %*% V[,h] 
    

    ### Step 2-3: Loading computation  ### 
    
    #Loading computation
    if(MODE=="PLS-SB"){
      # Loadings are set to be the directions
      P[,h]<-W[,h] 
      Q[,h]<-V[,h] 
    }else{
      tt = sum(T[,h]^2)
      uu = sum(U[,h]^2)
      P[,h] <- t(EE) %*% T[,h] / (tt)
      Q[,h] <- t(FF) %*% U[,h] / (uu)
    }
    
    ### Step 4: Regression (only in PLS_12) + Deflation    ### 
    
    #X-deflation
    EE <- EE - T[,h] %*% t(P[,h])
    
    #Y-deflation
    if(MODE=="PLS-R"){
      B[h,h] <- t(U[,h]) %*% T[,h] / (tt) 
      ## FF <- FF - B[h,h] * T[,h] %*% t(V[,h])
    }else{
      FF <- FF - U[,h]%*% t(Q[,h])  
    }
    
    # Saving intermediate results
    XX[[h+1]] = EE
    YY[[h+1]] = FF
    
    #X estimate
    X_hat[[h]] <- T[,1:h] %*% t(P[,1:h]) + X_MEAN
    
    #Y (and Beta) estimates
    if(MODE == "PLS-R"){
      Beta[[h]] <- W[,1:h] %*% solve(t(P[,1:h]) %*% W[,1:h] , B[1:h,1:h] %*% t(V[,1:h]))
      Y_hat[[h]] <- X_c %*% Beta[[h]] + Y_MEAN
    }else{
      Y_hat[[h]] <- U[,1:h] %*% t(Q[,1:h]) + Y_MEAN
    }
    
  }
  
  return(list(X_space_directions_locs = W,
              Y_space_directions = V,
              X_latent_scores = T,
              Y_latent_scores = U,
              X_loadings_locs = P,
              Y_loadings = Q,
              B = B,
              Beta_locs = Beta,
              X_mean_locs = X_mean,
              Y_mean = Y_mean,
              Y_hat = Y_hat,
              X_hat_locs = X_hat))
}
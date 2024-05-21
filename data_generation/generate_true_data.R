# % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
# % % Test: True data generation % %
# % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

rm(list = ls())
graphics.off()


## libraries ----

## algebraic utils
suppressMessages(library(pracma))

## statistical utils
suppressMessages(library(MASS))

## visualization
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))

## sampling
suppressMessages(library(sf))
suppressMessages(library(sp))
suppressMessages(library(raster))

## fdaPDE
suppressMessages(library(fdaPDE))
suppressMessages(library(fdaPDE2))


## sources ----
source("src/utils/directories.R")
source("src/utils/domain_and_locations.R")
source("src/utils/meshes.R")
source("src/utils/plots.R")
source("src/utils/errors.R")
source("src/utils/pls.R")


## functions ----
comp_sin <- function(locs, i){
  k1 <- c(1, 1, 2, 2, 1, 3, 2, 3, 1, 4, 2, 4, 3, 4, 4)
  k2 <- c(1, 2, 1, 2, 3, 2, 3, 3, 4, 1, 4, 2, 4, 3, 4)
  phi1 <- c(0, 0.6, -0.3, 0.2)
  phi2 <- c(0.35, 0, 0.2, 0.45)
  res<- sin(k1[i]*pi*(locs[,1]+phi1[i]))*sin(k2[i]*pi*(locs[,2]+phi2[i]))
  return(res)
}


## global options ----

path_data <- paste("data/fPLS/")
mkdir(path_data)


## data options ----
mode <- "PLS-R"
name_mesh <- "unit_square"
Beta_index <- 5
n_comp <- 4
sd_scores <- c(1, 0.8, 0.6, 0.5)


## generation accuracy parameters ----
n_nodes <- 50^2
n_stat_units <- 1000
seed <- 0


## domain & locations ----

## load domain
generated_domain <- generate_domain(name_mesh, n_nodes)
domain <- generated_domain$domain
domain_boundary <- generated_domain$domain_boundary
mesh <- generated_domain$mesh

## nodes
nodes <- mesh$nodes

## dimensions
n_nodes <- nrow(nodes)


## X loadings ----

## generating the loadings functions
X_loadings_true <- matrix(0, nrow = n_nodes, ncol = n_comp)
for (i in 1:n_comp) {
  X_loadings_true[, i] <- comp_sin(nodes, i)
  norm <- 1 # RMSE(X_loadings_true[, i])
  X_loadings_true[, i] <- X_loadings_true[, i] / norm
}

## computing the scores sd. (sigma)
## sd = sqrt(Var[score*loading]) = sqrt(Var[score] * semi_range(loading)^2) = sqrt(Var[score]) * semi_range(loading)
## sigma = sqrt(Var[score])
## => sigma = sd / semi_range(loading)
sd_s <- c(1, 0.8, 0.3, 0.2)
semi_range_X_loadings_true <- 0.5 * (apply(X_loadings_true, 2, max) - apply(X_loadings_true, 2, min))
sigma_s <- sd_s / semi_range_X_loadings_true


## Y loadings ----

if(Beta_index == 1)
  Y_loadings_true <- matrix(c(1, 0, 0, 0,
                              0, 1, 0, 0),
                            ncol = 4, byrow = TRUE)
if(Beta_index == 2)
  Y_loadings_true <- matrix(c(0, 1, 0, 0,
                              0, 0, 1, 0),
                            ncol = 4, byrow = TRUE)
if(Beta_index == 3)
  Y_loadings_true <- matrix(c(0, 0, 1, 0,
                              0, 0, 0, 1),
                            ncol = 4, byrow = TRUE)
if(Beta_index == 4)
  Y_loadings_true <- matrix(c(0, 0, 0, 1,
                              1, 0, 0, 0),
                            ncol = 4, byrow = TRUE)
if(Beta_index == 5)
  Y_loadings_true <- matrix(c(1, -3, 0, 8,
                              0, 3, 8, 0),
                            ncol = 4, byrow = TRUE)
if(Beta_index == 6)
  Y_loadings_true <- matrix(c(0, 0, 8, 0,
                              1, -3, 0, 5),
                            ncol = 4, byrow = TRUE)
n_resp <- nrow(Y_loadings_true)


## Sampling ----

set.seed(seed)
Sampled <- mvrnorm(10000, mu = rep(0, n_comp), diag(c(sigma_s^2)))
scores_true <- scale(Sampled, scale = FALSE)


## Data ----

X_true <- scores_true %*% t(X_loadings_true)
Y_true <- scores_true %*% t(Y_loadings_true)


## PLS ----

expected_results <- PLS(Y_true, X_true, n_comp = n_comp, mode = mode)


## save results ----

grid <- nodes
X_space_directions_true_grid <- expected_results$X_space_directions
Y_space_directions_true <- expected_results$Y_space_directions
X_loadings_true_grid <- expected_results$X_loadings
Y_loadings_true <- expected_results$Y_loadings
Beta_true_grid <- expected_results$Beta

sd_X_latent_scores <- sqrt(diag(cov(expected_results$X_latent_scores)))
sd_X_latent_scores/sd_X_latent_scores[1]

grid <- mesh

expected_results$Y_loadings

save(## data options
     mode,
     name_mesh,
     Beta_index,
     n_comp,
     ## generated data
     grid,
     X_space_directions_true_grid,
     Y_space_directions_true,
     X_loadings_true_grid,
     Y_loadings_true,
     sd_X_latent_scores,
     Beta_true_grid,
     ## file
     file = paste(path_data, paste(mode, name_mesh, "b", Beta_index, "comp", n_comp, sep = "_"), ".RData", sep = "")
)
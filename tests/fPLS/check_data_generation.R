# % %%%%%%%%%%%%%% %
# % % Test: fPLS % %
# % %%%%%%%%%%%%%% %

rm(list = ls())
graphics.off()

## global variables ----

test_suite <- "fPLS"
TEST_SUITE <- "fPLS"


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


## sources ----
source("src/utils/domain_and_locations.R")
source("src/utils/meshes.R")
source("src/utils/plots.R")
source("src/utils/wrappers.R")
source("src/utils/errors.R")
source("src/utils/pls.R")
source(paste("tests/", test_suite, "/utils/generate_2D_data.R", sep = ""))


## options ----

name_mesh <- "unit_square"
n_nodes <- 30^2
n_stat_units <- 50
locs_eq_nodes <- TRUE
n_comp <- 4
NSR_X_last_comp <- 0.1
NSR_Y <- 0.5
Beta_index <- 5
pls_mode <- "PLS-R"


## domain & locations ----

## load domain
generated_domain <- generate_domain(name_mesh, n_nodes)
domain <- generated_domain$domain
domain_boundary <- generated_domain$domain_boundary
loadings_true_generator <- generated_domain$loadings_true_generator
mesh <- generated_domain$mesh

## locations
locations <- generate_locations(generated_domain, locs_eq_nodes)


## data ----

## generate data
generated_data <- generate_2D_data(
  ## true data identificator
  pls_mode,     
  name_mesh,
  Beta_index,      
  n_comp,
  ## simulation setting
  domain,
  locs = locations,
  n_stat_units = n_stat_units,
  n_nodes = n_nodes,
  NSR_X_last_comp = NSR_X_last_comp,
  NSR_Y = NSR_Y,
  seed = 0,
  MODE = "PLS-R",
  VERBOSE = TRUE
)

## view X space directions
title <- "X space directions"
labels_cols <- as.list(paste("w", 1:n_comp, sep = ""))
X_space_directions_true_locs <- generated_data$expected_results$X_space_directions_true_locs
plot_list <- list()
for (h in 1:n_comp) {
  plot_list[[h]] <- plot.field_points(
    locations, X_space_directions_true_locs[, h],
    LEGEND = FALSE,
    boundary = domain_boundary, size = 1.5
  ) +
    standard_plot_settings_fields()
}
plot <- arrangeGrob(grobs = plot_list, nrow = 1)
plot <- labled_plots_grid(
  plot,
  title = title,
  labels_cols = labels_cols,
  width = 5,
  height = 5
)
grid.arrange(plot)

## view loadings
title <- "X loadings"
labels_cols <- as.list(paste("f", 1:n_comp, sep = ""))
X_loadings_true_locs <- generated_data$expected_results$X_loadings_true_locs
plot_list <- list()
for (h in 1:n_comp) {
  plot_list[[h]] <- plot.field_points(
    locations, X_loadings_true_locs[, h],
    LEGEND = FALSE,
    boundary = domain_boundary, size = 1.5
  ) +
    standard_plot_settings_fields()
}
plot <- arrangeGrob(grobs = plot_list, nrow = 1)
plot <- labled_plots_grid(
  plot,
  title = title,
  labels_cols = labels_cols,
  width = 5,
  height = 5
)
grid.arrange(plot)

## view scores
X_latent_scores_true <- generated_data$expected_results$X_latent_scores_true
colnames(X_latent_scores_true) <- c("t1", "t2", "t3", "t4")
boxplot(X_latent_scores_true, main = "Scores distributions")
grid()
boxplot(X_latent_scores_true, add = TRUE)
points(1:4, generated_data$sd_X_latent_scores, type = "l", lty = 2, col = "red")

pairs(X_latent_scores_true)

## noise distribution
noise <- generated_data$data$X - generated_data$expected_results$X_true_locs
boxplot(noise, xaxt = "n", main = "Noise distribution")
abline(h = 0, col = "red", lty = 1)
abline(h = c(-generated_data$sigma_noise_x, generated_data$sigma_noise_x), col = "red", lty = 2)

## noise spatial correlation
plot_list <- list()
for (i in 1:9) {
  plot_list[[i]] <- plot.field_points(locations, noise[i, ], boundary = domain_boundary) + standard_plot_settings_fields() + ggtitle(paste("Noise sample", i))
}
plot <- arrangeGrob(grobs = plot_list, nrow = 3)
grid.arrange(textGrob("Noise samples", gp = gpar(fontsize = 14, fontface = "bold")), plot, heights = c(1, 4 * 3))

## data samples
n_samples <- 3
title <- "X samples"
labels_cols <- list("True", "True at locations", "Observed")
labels_rows <- as.list(paste("x", 1:n_samples, sep = ""))
plot_list <- list()
for (i in 1:n_samples) {
  limits <- range(generated_data$data$X[i, ])
  plot_list[[(i - 1) * n_samples + 1]] <- plot.field_tile(domain$nodes, generated_data$expected_results$X_true[i, ], boundary = domain_boundary, limits = limits, ISOLINES = TRUE) + standard_plot_settings_fields()
  plot_list[[(i - 1) * n_samples + 2]] <- plot.field_points(locations, generated_data$expected_results$X_true_locs[i, ], boundary = domain_boundary, limits = limits, size = 1.5) + standard_plot_settings_fields()
  plot_list[[(i - 1) * n_samples + 3]] <- plot.field_points(locations, generated_data$data$X[i, ], boundary = domain_boundary, limits = limits, size = 1.5) + standard_plot_settings_fields()
}
plot <- arrangeGrob(grobs = plot_list, nrow = n_samples)
plot <- labled_plots_grid(plot, title, labels_cols, labels_rows, 5, 5)
grid.arrange(plot)

## centered data samples
n_samples <- 3
title <- "X centered samples"
labels_cols <- list("True", "True at locations", "Observed")
labels_rows <- as.list(paste("x", 1:n_samples, sep = ""))
plot_list <- list()
for (i in 1:n_samples) {
  limits <- range(generated_data$data$X[i, ] - generated_data$expected_results$X_mean_true_locs)
  plot_list[[(i - 1) * n_samples + 1]] <- plot.field_tile(domain$nodes, generated_data$expected_results$X_c_true[i, ], boundary = domain_boundary, limits = limits, ISOLINES = TRUE) + standard_plot_settings_fields()
  plot_list[[(i - 1) * n_samples + 2]] <- plot.field_points(locations, generated_data$expected_results$X_c_true_locs[i, ], boundary = domain_boundary, limits = limits, size = 1.5) + standard_plot_settings_fields()
  plot_list[[(i - 1) * n_samples + 3]] <- plot.field_points(locations, generated_data$data$X[i, ] - generated_data$expected_results$X_mean_true_locs, boundary = domain_boundary, limits = limits, size = 1.5) + standard_plot_settings_fields()
}
plot <- arrangeGrob(grobs = plot_list, nrow = n_samples)
plot <- labled_plots_grid(plot, title, labels_cols, labels_rows, 5, 5)
grid.arrange(plot)


## X last component
n_samples <- 3
title <- "X last component"
labels_cols <- list("True", "True at locations", "Observed")
labels_rows <- as.list(paste("x", 1:n_samples, sep = ""))
plot_list <- list()
for (i in 1:n_samples) {
  limits <- range(generated_data$expected_results$X_latent_scores_true[i, 4] *
                    generated_data$expected_results$X_loadings_true_locs[, 4] +
                    noise[i,])
  plot_list[[(i - 1) * n_samples + 1]] <- plot.field_tile(domain$nodes,
                                                          generated_data$expected_results$X_latent_scores_true[i, 4] *
                                                            generated_data$expected_results$X_loadings_true[, 4],
                                                          boundary = domain_boundary, limits = limits, ISOLINES = TRUE) + standard_plot_settings_fields()
  plot_list[[(i - 1) * n_samples + 2]] <- plot.field_points(locations,
                                                            generated_data$expected_results$X_latent_scores_true[i, 4] *
                                                              generated_data$expected_results$X_loadings_true_locs[, 4],
                                                            boundary = domain_boundary, limits = limits, size = 1.5) + standard_plot_settings_fields()
  plot_list[[(i - 1) * n_samples + 3]] <- plot.field_points(locations,
                                                            generated_data$expected_results$X_latent_scores_true[i, 4] *
                                                              generated_data$expected_results$X_loadings_true_locs[, 4] +
                                                              noise[i,],
                                                            boundary = domain_boundary, limits = limits, size = 1.5) + standard_plot_settings_fields()
}
plot <- arrangeGrob(grobs = plot_list, nrow = n_samples)
plot <- labled_plots_grid(plot, title, labels_cols, labels_rows, 5, 5)
grid.arrange(plot)

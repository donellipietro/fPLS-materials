# % %%%%%%%%%%%%%%%%%%%%%% %
# % % Test: TEST SUITE 1 % %
# % %%%%%%%%%%%%%%%%%%%%%% %

rm(list = ls())
graphics.off()

# setwd("")

## global variables ----

test_suite <- "test_suite_1"
TEST_SUITE <- "Test Suite 1"


# libraries ----

## fda
suppressMessages(library(fdaPDE))
suppressMessages(library(fdaPDE2))

## algebraic utils
suppressMessages(library(pracma))

## statistical utils
suppressMessages(library(MASS))

## data visualization
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))

## json
suppressMessages(library(jsonlite))

## sampling
suppressMessages(library(sf))
suppressMessages(library(sp))
suppressMessages(library(raster))


# sources ----

source("src/utils/cat.R")
source("src/utils/directories.R")
source("src/utils/meshes.R")
source("src/utils/domain_and_locations.R")
source("src/utils/wrappers.R")
source("src/generate_2D_data.R")
source("src/utils/errors.R")
source("src/utils/results_management.R")
source("src/utils/plots.R")

source(paste("tests/", test_suite, "/utils/models_evaluation.R", sep = ""))


# paths ----

path_results <- paste("results/", test_suite, "/", sep = "")
path_images <- paste("images/", test_suite, "/", sep = "")
mkdir(c(path_results, path_images))

path_queue <- paste("queue/", sep = "")


# files ----

file_log <- "log.txt"


# options ----

## execution flow modifiers
RUN <- list()
RUN$tests <- TRUE
RUN$plots <- TRUE

## global variables
RSTUDIO <- FALSE


## calibration parameters ----

seed <- 100 # for gcv
lambda_grid <- fdaPDE2::hyperparameters(10^seq(-6, 2, by = 0.2))


## visualization options ----

## colors used in the plots
colors <- c(brewer.pal(3, "Greens")[2], brewer.pal(3, "Blues")[2], brewer.pal(3, "Purples")[2])

## names and labels
names_models <- c("MV_PCA", "fPCA_off", "fPCA_gcv")
lables_models <- c("MV-PCA", "fPCA (no calibration)", "fPCA (gcv calibration)")

## resolution of the high resolution grid
n_nodes_HR_grid <- 1000


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Test: TEST SUITE 1 ----
cat.script_title(paste("Test:", TEST_SUITE))


## options ----
cat.section_title("Options")

## check arguments passed by terminal, set default if not present
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  RSTUDIO <- TRUE
  source(paste("tests/", test_suite, "/utils/generate_options.R", sep = ""))
  args[1] <- "test1"
  generate_options(args[1], path_queue)
  args[2] <- sort(list.files(path_queue))[1]
}

## read arguments provided
name_main_test <- args[1]
file_options <- args[2]

## load options
parsed_json <- fromJSON(paste(path_queue, file_options, sep = ""))
name_test <- parsed_json$test$name_test
name_mesh <- parsed_json$mesh$name_mesh
n_nodes <- parsed_json$dimensions$n_nodes
n_locs <- parsed_json$dimensions$n_locs
n_stat_units <- parsed_json$dimensions$n_stat_units
n_comp <- parsed_json$dimensions$n_comp
n_reps <- parsed_json$dimensions$n_reps
locs_eq_nodes <- parsed_json$data$locs_eq_nodes
NSR_X <- parsed_json$noise$NSR_X
# ... (add options if necessary)

## options visualization
cat.json(parsed_json)

## create results directory
path_results <- paste(path_results, name_main_test, "/", sep = "")
mkdir(path_results)
path_results <- paste(path_results, name_test, "/", sep = "")
mkdir(path_results)

## create images directory
path_images <- paste(path_images, name_main_test, "/", sep = "")
mkdir(path_images)
path_images <- paste(path_images, "single_tests/", sep = "")
mkdir(path_images)


## test ----
cat.section_title("Test")

## load domain
generated_domain <- generate_domain(name_mesh, n_nodes)
domain <- generated_domain$domain
boundary_domain <- generated_domain$boundary
loadings_true_generator <- generated_domain$loadings_true_generator
mesh <- generated_domain$mesh

## locations
locations <- generate_locations(generated_domain, locs_eq_nodes)

## fit the models n_reps times
if (RUN$tests) {
  for (i in 1:n_reps) {
    ## message
    cat(paste("\nBatch ", i, ":\n", sep = ""))

    ## generate data
    cat("- Data generation\n")
    generated_data <- generate_2D_data(
      domain = domain,
      locs = locations,
      n_stat_units = n_stat_units,
      n_comp = n_comp,
      seed = i,
      NSR_X = NSR_X,
      loadings_true_generator = loadings_true_generator,
    )

    ## assembly functional data
    data <- functional_data(
      domain = domain,
      locations = locations,
      X = generated_data$X
    )

    ## fit models
    source(paste("tests/", test_suite, "/templates/fit_and_evaluate.R", sep = ""))
  }
} else { ## end run tests
  cat("Skipped, relying on the saved results!\n")
}


## results analysis ----
cat.section_title("Results analysis")


## laod saved results
if (RUN$plots) {
  source(paste("tests/", test_suite, "/templates/load_results.R", sep = ""))
}


### quantitative analysis ----
cat.subsection_title("Quantitative analysis")

## open a pdf where to save plots (quantitative analysis)
if (!RSTUDIO) {
  pdf(file = paste(path_images, name_test, "_quantitative.pdf", sep = ""))
}

## plots
if (RUN$plots) {
  ## Time
  plot <- plot.grouped_boxplots(
    times,
    values_name = "Time [seconds]",
    group_name = "",
    group_labels = "",
    subgroup_name = "Approaches",
    subgroup_labels = lables_models,
    subgroup_colors = colors
  ) + standard_plot_settings() + ggtitle("Time")
  print(plot)

  ## RMSE
  names <- c("reconstruction", "reconstruction_locs")
  titles <- c("Reconstruction", "Reconstruction at locations")
  for (i in 1:length(names)) {
    name <- names[i]
    title <- titles[i]
    plot <- plot.grouped_boxplots(
      rmses[[name]],
      values_name = "RMSE",
      group_name = "",
      group_labels = "",
      subgroup_name = "Approaches",
      subgroup_labels = lables_models,
      subgroup_colors = colors
    ) + standard_plot_settings() + ggtitle(title)
    print(plot)
  }
  names <- c("loadings", "loadings_locs", "scores")
  titles <- c("Loadings", "Loadings at locations", "Scores")
  for (i in 1:length(names)) {
    name <- names[i]
    title <- titles[i]
    plot <- plot.grouped_boxplots(
      rmses[[name]],
      values_name = "RMSE",
      subgroup_name = "Approaches",
      subgroup_labels = lables_models,
      subgroup_colors = colors
    ) + standard_plot_settings() + ggtitle(title)
    print(plot)
  }
} ## end plot results

## close pdf quantitative
dev.off()


### qualitative analysis ----
cat.subsection_title("Qualitative analysis")

## open a pdf where to save plots (qualitative analysis)
if (!RSTUDIO) {
  pdf(file = paste(path_images, name_test, "_qualitative.pdf", sep = ""), width = 10, height = 7)
}

## plots
if (RUN$plots) {
  selected <- c(1, 2, 3)
  names_models_plot <- names_models[selected]
  lables_models_plot <- lables_models[selected]
  colors_plot <- colors[selected]
  labels_cols <- c("True", "Quantile 0", "Quantile 0.5", "Quantile  1")
  labels_rows <- c("f1", "f2", "f3")

  plots <- list()
  plots_locs <- list()
  plots_HR_clean <- list()
  plots_HR <- list()

  ## generate figures
  for (m in 1:length(names_models_plot)) {
    name_model <- names_models_plot[m]
    lable_model <- lables_models_plot[m]
    plot_list <- list()
    plot_list_locs <- list()
    plot_list_HR_clean <- list()
    plot_list_HR <- list()
    indexes <- tapply(rmses$loadings[[name_model]], rmses$loadings$Group, function(x) {
      match(quantile(x, c(0, 0.5, 1)), x)
    })
    for (i in 1:n_comp) {
      limits <- range(loadings_true_HR[, i])
      breaks <- seq(limits[1], limits[2], length = 10)
      plot_list[[4 * (i - 1) + 1]] <- plot.field_points(domain$nodes, loadings_true[, i], boundary = domain_boundary, size = 1.5, LEGEND = FALSE) +
        standard_plot_settings_fields()
      plot_list_locs[[4 * (i - 1) + 1]] <- plot.field_points(locations, loadings_true_locs[, i], boundary = domain_boundary, size = 1.5, LEGEND = FALSE) +
        standard_plot_settings_fields()
      plot_list_HR_clean[[4 * (i - 1) + 1]] <- plot.field_tile(grid, loadings_true_HR[, i], boundary = domain_boundary, LEGEND = FALSE) +
        standard_plot_settings_fields()
      plot_list_HR[[4 * (i - 1) + 1]] <- plot.field_tile(grid, loadings_true_HR[, i], boundary = domain_boundary, limits = limits, breaks = breaks, LEGEND = FALSE) +
        standard_plot_settings_fields()
      for (j in 1:3) {
        plot_list[[4 * (i - 1) + j + 1]] <- plot.field_points(domain$nodes, loadings[[name_model]][[indexes[[i]][j]]][, i], boundary = domain_boundary, size = 1.5) +
          standard_plot_settings_fields()
        plot_list_locs[[4 * (i - 1) + j + 1]] <- plot.field_points(locations, loadings_locs[[name_model]][[indexes[[i]][j]]][, i], boundary = domain_boundary, size = 1.5) +
          standard_plot_settings_fields()
        plot_list_HR_clean[[4 * (i - 1) + j + 1]] <- plot.field_tile(grid, loadings_HR[[name_model]][[indexes[[i]][j]]][, i], boundary = domain_boundary) +
          standard_plot_settings_fields()
        if (m != 3) {
          plot_list_HR[[4 * (i - 1) + j + 1]] <- plot.field_tile(grid, loadings_HR[[name_model]][[indexes[[i]][j]]][, i], boundary = domain_boundary, limits = limits, breaks = breaks) +
            standard_plot_settings_fields()
        } else {
          plot_list_HR[[4 * (i - 1) + j + 1]] <- plot.field_tile(grid, loadings_HR[[name_model]][[indexes[[i]][j]]][, i], boundary = domain_boundary, limits = limits) +
            standard_plot_settings_fields()
        }
      }
    }
    plots[[m]] <- arrangeGrob(grobs = plot_list, nrow = 3)
    plots[[m]] <- labled_plots_grid(plots[[m]], lable_model, labels_cols, labels_rows)
    plots_locs[[m]] <- arrangeGrob(grobs = plot_list_locs, nrow = 3)
    plots_locs[[m]] <- labled_plots_grid(plots_locs[[m]], lable_model, labels_cols, labels_rows)
    plots_HR_clean[[m]] <- arrangeGrob(grobs = plot_list_HR_clean, nrow = 3)
    plots_HR_clean[[m]] <- labled_plots_grid(plots_HR_clean[[m]], lable_model, labels_cols, labels_rows)
    plots_HR[[m]] <- arrangeGrob(grobs = plot_list_HR, nrow = 3)
    plots_HR[[m]] <- labled_plots_grid(plots_HR[[m]], lable_model, labels_cols, labels_rows)
  }

  for (m in 1:length(names_models_plot)) {
    grid.arrange(plots_locs[[m]])
  }
  for (m in 1:length(names_models_plot)) {
    grid.arrange(plots[[m]])
  }
  for (m in 1:length(names_models_plot)) {
    grid.arrange(plots_HR_clean[[m]])
  }
  for (m in 1:length(names_models_plot)) {
    grid.arrange(plots_HR[[m]])
  }
}

dev.off()


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## remove options file
file.remove(paste(path_queue, file_options, sep = ""))

cat("\n\n")

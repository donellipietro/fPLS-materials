# % %%%%%%%%%%%%%% %
# % % Test: fPLS % %
# % %%%%%%%%%%%%%% %

rm(list = ls())
graphics.off()


## global variables ----

test_suite <- "fPLS"
TEST_SUITE <- "fPLS"

mode_MV <- "PLS-A"
mode_fun <- "fPLS-A"


## prerequisite ----
source(paste("tests/", test_suite, "/utils/generate_options.R", sep = ""))


## libraries ----

## json
suppressMessages(library(jsonlite))

## algebraic utils
suppressMessages(library(pracma))

## data visualization
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))

## sources ----
source("src/utils/directories.R")
source("src/utils/results_management.R")
source("src/utils/plots.R")


## paths ----
path_options <- paste("queue/", sep = "")
path_results <- paste("results/", test_suite, "/", sep = "")
path_images <- paste("images/", test_suite, "/", sep = "")
file_log <- "log.txt"


## options ----

## colors used in the plots
colors <- c(brewer.pal(3, "Greys")[3], brewer.pal(3, "Blues")[2:3])

## names and labels
names_models <- c("MV_PLS", "fPLS_off", "fPLS_gcv")
lables_models <- c("MV-PLS", "fPLS (no calibration)", "fPLS (gcv calibration)")


## load data ----

## check arguments passed by terminal
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  args[1] <- "test1"
}

## main test name
name_main_test <- args[1]
cat(paste("\nTest selected:", name_main_test, "\n"))
path_results <- paste(path_results, mode_MV, "/", sep = "")
path_results <- paste(path_results, name_main_test, "/", sep = "")
path_images <- paste(path_images, mode_MV, "/", sep = "")
path_images <- paste(path_images, name_main_test, "/", sep = "")
mkdir(path_images)

## generate options
generate_options(name_main_test, path_options)

## list of available tests
file_test_vect <- sort(list.files(path_options))

## room for solutions
names_columns <- c("Group", "n_stat_units", "NSR_X_last_comp", "NSR_Y", names_models)
empty_df <- data.frame(matrix(NaN, nrow = 0, ncol = length(names_columns)))
colnames(empty_df) <- names_columns
times <- empty_df
rmses <- list()
irmses <- list()
angles <- list()


for (file_test in file_test_vect) {
  ## load specs
  file_json <- paste(path_options, file_test, sep = "")
  parsed_json <- fromJSON(file_json)
  name_test <- parsed_json$test$name_test
  n_stat_units <- parsed_json$dimensions$n_stat_units
  NSR_X_last_comp <- parsed_json$noise$NSR_X_last_comp
  NSR_Y <- parsed_json$noise$NSR_Y
  n_reps <- parsed_json$dimensions$n_reps
  
  cat(paste("\nTest ", name_test, ":\n", sep = ""))
  
  ## load batches
  for (i in 1:n_reps) {
    ## laod batch and log if not present
    sink(file_log, append = TRUE)
    tryCatch(
      {
        path_batch <- paste(path_results, name_test, "/", "batch_", i, "/", sep = "")
        load(paste(path_batch, "batch_", i, "_results_evaluation.RData", sep = ""))
      },
      error = function(e) {
        cat(paste("Error in test ", name_test, " - batch ", i, ": ", conditionMessage(e), "\n", sep = ""))
      }
    )
    sink()
    
    ## times
    times <- add_results(
      times,
      c(
        list(n_stat_units = n_stat_units, NSR_X_last_comp = NSR_X_last_comp, NSR_Y = NSR_Y),
        extract_new_results(results_evaluation, names_models, "execution_time")
      ),
      names_columns
    )
    
    ## rmse
    names <- c("Y_reconstruction", "X_reconstruction_locs")
    if(mode_MV == "PLS-R") {
      names <- c(names, "Beta_locs")
    }
    for (name in names) {
      rmses[[name]] <- add_results(
        rmses[[name]],
        c(
          list(n_stat_units = n_stat_units, NSR_X_last_comp = NSR_X_last_comp, NSR_Y = NSR_Y),
          extract_new_results(results_evaluation, names_models, c("rmse", name))
        ),
        names_columns
      )
    }
    
    names <- c("Y_space_directions", "X_space_directions_locs",
               "Y_latent_scores", "X_latent_scores",
               "Y_loadings", "X_loadings_locs")
    for (name in names) {
      rmses[[name]] <- add_results(
        rmses[[name]],
        c(
          list(n_stat_units = n_stat_units, NSR_X_last_comp = NSR_X_last_comp, NSR_Y = NSR_Y),
          extract_new_results(results_evaluation, names_models, c("rmse", name))
        ),
        names_columns
      )
    }
    
    cat(paste("- Batch", i, "loaded\n"))
  }
  
  ## remove option file
  file.remove(file_json)
}


## analysis ----

## names
name_aggregation_option_vect <- c("n_stat_units", "NSR_X_last_comp", "NSR_Y")
name_group_vect <- c("N", "NSR_X_lc", "NSR_Y")

### time complexity ----

## open a pdf where to save the plots
pdf(paste(path_images, "time_complexity.pdf", sep = ""), width = 20, height = 20)

## data and titles
data_plot <- times
values_name <- "Time [seconds]"
title_vect <- paste(
  "Execution times w.r.t the",
  c(
    "number of statistical units (N)",
    "NSR X last component", "NSR Y last component"
  )
)

## options
options_grid <- list()
for (name_ao in name_aggregation_option_vect) {
  options_grid[[name_ao]] <- unique(data_plot[, name_ao])
}


for (i in c(1)) {
  boxplot_list <- list()
  plot_list <- list()
  plot_loglog_list <- list()
  plot_loglog_normalized_list <- list()
  
  name_aggregation_option <- name_aggregation_option_vect[i]
  group_name <- name_group_vect[i]
  title <- title_vect[i]
  
  options_grid_selected <- options_grid
  options_grid_selected[[name_aggregation_option]] <- NULL
  names_options_selected <- names(options_grid_selected)
  labels_options_selected <- name_group_vect[-i]
  mg <- meshgrid(options_grid_selected[[1]], options_grid_selected[[2]])
  combinations_options <- data.frame(as.vector(mg[[1]]), as.vector(mg[[2]]))
  colnames(combinations_options) <- names_options_selected
  
  limits <- c(0, max(data_plot[, names_models]))
  range <- range(times[, names_models])
  
  labels_rows <- paste(labels_options_selected[1], "=", options_grid_selected[[1]])
  labels_cols <- paste(labels_options_selected[2], "=", options_grid_selected[[2]])
  
  for (j in 1:nrow(combinations_options)) {
    ## data preparation
    data_plot_trimmed <- data_plot[
      data_plot[, names_options_selected[1]] == combinations_options[j, 1] &
        data_plot[, names_options_selected[2]] == combinations_options[j, 2],
      c(name_aggregation_option, names_models)
    ]
    colnames(data_plot_trimmed) <- c("Group", names_models)
    indexes <- which(!is.nan(colSums(data_plot_trimmed[, names_models])))
    data_plot_aggregated <- aggregate(. ~ Group, data = data_plot_trimmed, FUN = median)
    
    ## plots
    boxplot_list[[j]] <- plot.grouped_boxplots(
      data_plot_trimmed[, c("Group", names(indexes))],
      values_name = NULL, # values_name,
      group_name = group_name,
      subgroup_name = "Approches",
      subgroup_labels = lables_models[indexes],
      subgroup_colors = colors[indexes],
      limits = limits,
      LEGEND = FALSE
    ) + standard_plot_settings()
    plot_list[[j]] <- plot.multiple_lines(
      data_plot_aggregated[, c("Group", names(indexes))],
      values_name = NULL, # values_name,
      x_name = group_name,
      x_breaks = TRUE,
      subgroup_name = "Approches",
      subgroup_labels = lables_models[indexes],
      subgroup_colors = colors[indexes],
      LEGEND = FALSE,
      limits = limits,
      NORMALIZED = FALSE,
      LOGX = TRUE
    ) + standard_plot_settings()
    plot_loglog_list[[j]] <- plot.multiple_lines(
      data_plot_aggregated[, c("Group", names(indexes))],
      values_name = NULL,
      x_name = group_name,
      x_breaks = TRUE,
      subgroup_name = "Approches",
      subgroup_labels = lables_models[indexes],
      subgroup_colors = colors[indexes],
      LEGEND = FALSE,
      limits = range,
      NORMALIZED = FALSE,
      LOGLOG = TRUE
    ) + standard_plot_settings()
    plot_loglog_normalized_list[[j]] <- plot.multiple_lines(
      data_plot_aggregated[, c("Group", names(indexes))],
      values_name = NULL,
      x_name = group_name,
      x_breaks = TRUE,
      subgroup_name = "Approches",
      subgroup_labels = lables_models[indexes],
      subgroup_colors = colors[indexes],
      LEGEND = FALSE,
      limits = c(1, limits[2]),
      NORMALIZED = TRUE,
      LOGLOG = TRUE
    ) + standard_plot_settings()
  }
  boxplot <- arrangeGrob(grobs = boxplot_list, ncol = length(labels_cols))
  boxplot <- labled_plots_grid(boxplot, title, labels_cols, labels_rows, 9, 7)
  grid.arrange(boxplot)
  plot <- arrangeGrob(grobs = plot_list, ncol = length(labels_cols))
  plot <- labled_plots_grid(plot, title, labels_cols, labels_rows, 9, 7)
  grid.arrange(plot)
  plot_loglog <- arrangeGrob(grobs = plot_loglog_list, ncol = length(labels_cols))
  plot_loglog <- labled_plots_grid(plot_loglog, title, labels_cols, labels_rows, 9, 7)
  grid.arrange(plot_loglog)
  plot_loglog_normalized <- arrangeGrob(grobs = plot_loglog_normalized_list, ncol = length(labels_cols))
  plot_loglog_normalized <- labled_plots_grid(plot_loglog_normalized, title, labels_cols, labels_rows, 9, 7)
  grid.arrange(plot_loglog_normalized)
}

dev.off()


### overall quantitative results ----


## open a pdf where to save the plots
pdf(paste(path_images, "overall_quantitative_results.pdf", sep = ""), width = 20, height = 20)

## Y reconstruction RMSE at locations
data_plot <- rmses[["Y_reconstruction"]][rmses[["Y_reconstruction"]]$Group == "4", ]
title_vect <- paste(
  "Y reconstruction RMSE at locations w.r.t the",
  c(
    "number of statistical units (N)",
    "NSR X last component", "NSR Y last component"
  )
)
limits <- NULL
source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))

## X reconstruction RMSE at locations
data_plot <- rmses[["X_reconstruction_locs"]][rmses[["X_reconstruction_locs"]]$Group == "4", ]
title_vect <- paste(
  "X reconstruction RMSE at locations w.r.t the",
  c(
    "number of statistical units (N)",
    "NSR X last component", "NSR Y last component"
  )
)
limits <- NULL
source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))

if(mode_MV == "PLS-R") {
  ## Beta RMSE at locations
  data_plot <- rmses[["Beta_locs"]][rmses[["Beta_locs"]]$Group == "4", ]
  title_vect <- paste(
    "Beta RMSE at locations w.r.t the",
    c(
      "number of statistical units (N)",
      "NSR X last component", "NSR Y last component"
    )
  )
  limits <- NULL
  source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))
}

dev.off()

### X quantitative results ----

## open a pdf where to save the plots
pdf(paste(path_images, "X_quantitative_results.pdf", sep = ""), width = 20, height = 20)

name_vect <- c("1", "2", "3", "4")
for (k in 1:4) {
  data_plot <- rmses[["X_space_directions_locs"]][rmses[["X_space_directions_locs"]]$Group == k, ]
  title_vect <- paste(
    "X space directions at locations RMSE component", name_vect[k], "w.r.t the",
    c(
      "number of nodes (K)",
      "number of statistical units (N)",
      "number of locations (S)"
    )
  )
  limits <- NULL
  source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))
}
for (k in 1:4) {
  data_plot <- rmses[["X_latent_scores"]][rmses[["X_latent_scores"]]$Group == k, ]
  title_vect <- paste(
    "X latent scores RMSE component", name_vect[k], "w.r.t the",
    c(
      "number of nodes (K)",
      "number of statistical units (N)",
      "number of locations (S)"
    )
  )
  limits <- NULL
  source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))
}
for (k in 1:4) {
  data_plot <- rmses[["X_loadings_locs"]][rmses[["X_loadings_locs"]]$Group == k, ]
  title_vect <- paste(
    "X loadings at locations RMSE component", name_vect[k], "w.r.t the",
    c(
      "number of nodes (K)",
      "number of statistical units (N)",
      "number of locations (S)"
    )
  )
  limits <- NULL
  source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))
}

dev.off()


### Y quantitative results ----

## open a pdf where to save the plots
pdf(paste(path_images, "Y_quantitative_results.pdf", sep = ""), width = 20, height = 20)

name_vect <- c("1", "2", "3", "4")
for (k in 1:4) {
  data_plot <- rmses[["Y_space_directions"]][rmses[["Y_space_directions"]]$Group == k, ]
  title_vect <- paste(
    "Y space directions RMSE component", name_vect[k], "w.r.t the",
    c(
      "number of nodes (K)",
      "number of statistical units (N)",
      "number of locations (S)"
    )
  )
  limits <- NULL
  source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))
}
for (k in 1:4) {
  data_plot <- rmses[["Y_latent_scores"]][rmses[["Y_latent_scores"]]$Group == k, ]
  title_vect <- paste(
    "Y latent scores RMSE component", name_vect[k], "w.r.t the",
    c(
      "number of nodes (K)",
      "number of statistical units (N)",
      "number of locations (S)"
    )
  )
  limits <- NULL
  source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))
}
for (k in 1:4) {
  data_plot <- rmses[["Y_loadings"]][rmses[["Y_loadings"]]$Group == k, ]
  title_vect <- paste(
    "Y loadings at locations RMSE component", name_vect[k], "w.r.t the",
    c(
      "number of nodes (K)",
      "number of statistical units (N)",
      "number of locations (S)"
    )
  )
  limits <- NULL
  source(paste("tests/", test_suite, "/templates/plot_overall.R", sep = ""))
}

dev.off()


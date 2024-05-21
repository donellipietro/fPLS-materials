## labels plot grids
labels_cols <- c("True", "Quantile 0", "Quantile 0.5", "Quantile  1")

## room for plots
plots <- list()
plots_locs <- list()
plots_HR_clean <- list()
plots_HR <- list()

## generate figures
for (m in 1:length(names_models)) {
  name_model <- names_models[m]
  lable_model <- lables_models[m]
  plot_list <- list()
  plot_list_locs <- list()
  plot_list_HR_clean <- list()
  plot_list_HR <- list()
  indexes <- tapply(rmses$X_space_directions_locs[[name_model]], rmses$X_space_directions_locs$Group, function(x) {
    match(quantile(x, c(0, 0.5, 1)), x)
  })
  for (i in 1:n_comp) {
    limits <- range(X_space_directions_true_HR[, i])
    breaks <- seq(limits[1], limits[2], length = 10)
    plot_list[[4 * (i - 1) + 1]] <- plot.field_points(domain$nodes, X_space_directions_true[, i], boundary = domain_boundary, size = 1.5, LEGEND = FALSE) +
      standard_plot_settings_fields()
    plot_list_locs[[4 * (i - 1) + 1]] <- plot.field_points(locations, X_space_directions_true_locs[, i], boundary = domain_boundary, size = 1.5, LEGEND = FALSE) +
      standard_plot_settings_fields()
    plot_list_HR_clean[[4 * (i - 1) + 1]] <- plot.field_tile(grid_HR, X_space_directions_true_HR[, i], boundary = domain_boundary, LEGEND = FALSE) +
      standard_plot_settings_fields()
    plot_list_HR[[4 * (i - 1) + 1]] <- plot.field_tile(grid_HR, X_space_directions_true_HR[, i], boundary = domain_boundary, limits = limits, breaks = breaks, LEGEND = FALSE) +
      standard_plot_settings_fields()
    for (j in 1:3) {
      plot_list[[4 * (i - 1) + j + 1]] <- plot.field_points(domain$nodes, X_space_directions[[name_model]][[indexes[[i]][j]]][, i], boundary = domain_boundary, size = 1.5) +
        standard_plot_settings_fields()
      plot_list_locs[[4 * (i - 1) + j + 1]] <- plot.field_points(locations, X_space_directions_locs[[name_model]][[indexes[[i]][j]]][, i], boundary = domain_boundary, size = 1.5) +
        standard_plot_settings_fields()
      plot_list_HR_clean[[4 * (i - 1) + j + 1]] <- plot.field_tile(grid_HR, X_space_directions_HR[[name_model]][[indexes[[i]][j]]][, i], boundary = domain_boundary) +
        standard_plot_settings_fields()
      plot_list_HR[[4 * (i - 1) + j + 1]] <- plot.field_tile(grid_HR, X_space_directions_HR[[name_model]][[indexes[[i]][j]]][, i], boundary = domain_boundary, limits = limits, breaks = breaks) +
        standard_plot_settings_fields()
    }
  }
  plots[[m]] <- arrangeGrob(grobs = plot_list, nrow = 4)
  plots[[m]] <- labled_plots_grid(plots[[m]], lable_model, labels_cols, labels_rows)
  plots_locs[[m]] <- arrangeGrob(grobs = plot_list_locs, nrow = 4)
  plots_locs[[m]] <- labled_plots_grid(plots_locs[[m]], lable_model, labels_cols, labels_rows)
  plots_HR_clean[[m]] <- arrangeGrob(grobs = plot_list_HR_clean, nrow = 4)
  plots_HR_clean[[m]] <- labled_plots_grid(plots_HR_clean[[m]], lable_model, labels_cols, labels_rows)
  plots_HR[[m]] <- arrangeGrob(grobs = plot_list_HR, nrow = 4)
  plots_HR[[m]] <- labled_plots_grid(plots_HR[[m]], lable_model, labels_cols, labels_rows)
}

for (m in 1:length(names_models)) {
  grid.arrange(plots_locs[[m]])
}
for (m in 1:length(names_models)) {
  grid.arrange(plots[[m]])
}
for (m in 1:length(names_models)) {
  grid.arrange(plots_HR_clean[[m]])
}
for (m in 1:length(names_models)) {
  grid.arrange(plots_HR[[m]])
}

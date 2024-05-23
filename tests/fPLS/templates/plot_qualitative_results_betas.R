## labels plot grids
labels_cols <- c("True", "Quantile 0", "Quantile 0.5", "Quantile  1")

for(beta_idx in c(1,2)){
  
  idx=as.character(beta_idx)
  labels_rows <- c(paste0("beta",idx," - 1"),paste0("beta",idx," - 1:2"), paste0("beta",idx," - 1:3"), paste0("beta",idx," - 1:4"))
  
  ## room for plots
  plots <- list()
  plots_locs <- list()
  plots_HR_clean <- list()
  plots_HR <- list()
  
  plots_true <- list()
  plots_locs_true <- list()
  plots_HR_clean_true <- list()
  plots_HR_true <- list()
  
  n_rows=4
    
  ## generate figures
  for (m in 1:length(names_models)) {
    name_model <- names_models[m]
    lable_model <- lables_models[m]
    
    plot_list <- list()
    plot_list_locs <- list()
    plot_list_HR_clean <- list()
    plot_list_HR <- list()
    
    plot_list_true <- list()
    plot_list_locs_true <- list()
    plot_list_HR_clean_true <- list()
    plot_list_HR_true <- list()
    
    indexes <- tapply(rmses$Beta_locs[[name_model]], rmses$Beta_locs$Group, function(x) {
      match(quantile(x, c(0, 0.5, 1)), x)
    })
    
    
    limits <- range(Betas_true_HR[, beta_idx])
    breaks <- seq(limits[1], limits[2], length = 10)


    for (i in 1:n_comp) { ##cambiato qui
      

      limits <- range(Betas_true_HR[, beta_idx])
      
       plot_list[[n_rows * (i - 1) + 1]] <- plot.field_points(domain$nodes, Betas_true[, beta_idx], boundary = domain_boundary, size = 1.5, LEGEND = FALSE) +
         standard_plot_settings_fields()
       
       plot_list_locs[[n_rows * (i - 1) + 1]] <- plot.field_points(locations, Betas_true_locs[, beta_idx], boundary = domain_boundary, size = 1.5, LEGEND = FALSE) +
         standard_plot_settings_fields()
       
       plot_list_HR_clean[[n_rows * (i - 1) + 1]] <- plot.field_tile(grid_HR, Betas_true_HR[, beta_idx], boundary = domain_boundary, LEGEND = FALSE) +
         standard_plot_settings_fields()
       
       plot_list_HR[[n_rows * (i - 1) + 1]] <- plot.field_tile(grid_HR, Betas_true_HR[, beta_idx], boundary = domain_boundary, limits = limits, breaks = breaks, LEGEND = FALSE) +
         standard_plot_settings_fields()
       
        

        for (j in 1:3) {
        plot_list[[n_rows * (i - 1) + j + 1]] <- plot.field_points(domain$nodes, Betas[[name_model]][[indexes[[i]][j]]][[beta_idx]][, i], boundary = domain_boundary, size = 1.5) +
          standard_plot_settings_fields()
        plot_list_locs[[n_rows * (i - 1) + j + 1]] <- plot.field_points(locations, Betas_locs[[name_model]][[indexes[[i]][j]]][[beta_idx]][, i], boundary = domain_boundary, size = 1.5) +
          standard_plot_settings_fields()
        plot_list_HR_clean[[n_rows * (i - 1) + j + 1]]<- plot.field_tile(grid_HR, Betas_HR[[name_model]][[indexes[[i]][j]]][[beta_idx]][, i], boundary = domain_boundary) +
          standard_plot_settings_fields()
        plot_list_HR[[n_rows * (i - 1) + j + 1]] <- plot.field_tile(grid_HR, Betas_HR[[name_model]][[indexes[[i]][j]]][[beta_idx]][, i], boundary = domain_boundary, limits = limits, breaks = breaks) +
          standard_plot_settings_fields()
        }
        
      }
      
      plots[[m]] <- arrangeGrob(grobs = plot_list, nrow = n_rows)
      plots[[m]] <- labled_plots_grid(plots[[m]], lable_model, labels_cols, labels_rows)
      
      plots_locs[[m]] <- arrangeGrob(grobs = plot_list_locs, nrow = n_rows)
      plots_locs[[m]] <- labled_plots_grid(plots_locs[[m]], lable_model, labels_cols, labels_rows)
      
      plots_HR_clean[[m]] <- arrangeGrob(grobs = plot_list_HR_clean, nrow = n_rows)
      plots_HR_clean[[m]] <- labled_plots_grid(plots_HR_clean[[m]], lable_model, labels_cols, labels_rows)
      
      plots_HR[[m]] <- arrangeGrob(grobs = plot_list_HR, nrow = n_rows)
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

}

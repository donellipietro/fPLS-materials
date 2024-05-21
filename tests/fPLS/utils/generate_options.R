generate_options <- function(name_main_test, path_queue) {
  ## create the directory if it does not exist
  mkdir(c(path_queue))
  
  ## name test suite
  test_suite <- "fPLS"
  
  ## generate the options json files
  switch(name_main_test,
         test1 = {
           
           # Test 1:
           
           ## options grid
           n_stat_units_vect <- c(50, 100, 200)
           NSR_X_last_comp_vect <- seq(0.1, 0.3, by = 0.2)
           
           ## fixed options
           name_mesh <- "unit_square"
           name_mesh_short <- "us"
           n_nodes <- 30^2
           locs_eq_nodes <- TRUE
           n_comp <- 4
           Beta_index <- 5
           NSR_Y <- 0.1
           n_reps <- 3
           mean <- FALSE
           
           ## assembly jsons           
           
           for (n_stat_units in n_stat_units_vect) {
             for (NSR_X_last_comp in NSR_X_last_comp_vect) {
               name_test <- paste(
                 test_suite, name_main_test, name_mesh_short,
                 sprintf("%0*d", 4, n_nodes),
                 sprintf("%0*d", 4, n_nodes), ## locations = nodes
                 sprintf("%0*d", 4, n_stat_units),
                 Beta_index,
                 NSR_X_last_comp, NSR_Y, n_comp,
                 sep = "_"
               )
               json_list <- list(
                 test = list(
                   name_test = name_test
                 ),
                 mesh = list(
                   name_mesh = name_mesh
                 ),
                 dimensions = list(
                   n_nodes = n_nodes,
                   n_locs = n_nodes, ## locations = nodes
                   n_stat_units = n_stat_units,
                   n_comp = n_comp,
                   n_reps = n_reps
                 ),
                 data = list(
                   mean = mean,
                   locs_eq_nodes = locs_eq_nodes,
                   Beta_index = Beta_index
                 ),
                 noise = list(
                   NSR_X_last_comp = NSR_X_last_comp,
                   NSR_Y = NSR_Y
                 )
               )
               write_json(json_list, path = paste(path_queue, name_test, ".json", sep = ""))
             }
           }
         },
         test2 = {
           
           # Test 2:
           
           ## options grid
           n_locs_vect <- c(400, 600, 800, 1000)
           n_stat_units_vect <- c(50, 100, 200)
           NSR_X_last_comp_vect <- seq(0.1, 1.1, by = 0.2)
           
           ## fixed options
           name_mesh <- "unit_square"
           name_mesh_short <- "us"
           n_nodes <- 30^2
           locs_eq_nodes <- FALSE
           n_comp <- 4
           Beta_index <- 5
           NSR_Y <- 0.5
           n_reps <- 3
           mean <- FALSE
           
           ## assembly jsons
           for (n_locs in n_locs_vect) {
             for (n_stat_units in n_stat_units_vect) {
               for (NSR_X_last_comp in NSR_X_last_comp_vect) {
                 name_test <- paste(
                   test_suite, name_main_test, name_mesh_short,
                   sprintf("%0*d", 4, n_nodes),
                   sprintf("%0*d", 4, n_locs),
                   sprintf("%0*d", 4, n_stat_units),
                   Beta_index,
                   NSR_X_last_comp, NSR_Y, n_comp,
                   sep = "_"
                 )
                 json_list <- list(
                   test = list(
                     name_test = name_test
                   ),
                   mesh = list(
                     name_mesh = name_mesh
                   ),
                   dimensions = list(
                     n_nodes = n_nodes,
                     n_locs = n_locs,
                     n_stat_units = n_stat_units,
                     n_comp = n_comp,
                     n_reps = n_reps
                   ),
                   data = list(
                     mean = mean,
                     locs_eq_nodes = locs_eq_nodes,
                     Beta_index = Beta_index
                   ),
                   noise = list(
                     NSR_X_last_comp = NSR_X_last_comp,
                     NSR_Y = NSR_Y
                   )
                 )
                 write_json(json_list, path = paste(path_queue, name_test, ".json", sep = ""))
               }
             }
           }
         },
         {
           stop(paste("The test", name_main_test, "does not exist"))
         }
  )
}

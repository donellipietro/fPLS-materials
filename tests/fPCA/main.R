# % %%%%%%%%%%%%%% %
# % % Test: fPCA % %
# % %%%%%%%%%%%%%% %

rm(list = ls())
graphics.off()

# setwd("")

## global variables ----

test_suite <- "fPCA"
TEST_SUITE <- "fPCA"


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

## force testing even if a fit is already available
FORCE <- TRUE

## execution flow modifiers
RUN <- list()
RUN$tests <- TRUE
RUN$analysis <- TRUE
RUN$quantitative_analysis <- TRUE
RUN$qualitative_analysis <- TRUE

## global variables
RSTUDIO <- FALSE


## calibration parameters ----

seed <- 100 # for gcv
lambda_grid <- fdaPDE2::hyperparameters(10^seq(-6, 2, by = 0.2))


## visualization options ----

## colors used in the plots
colors <- c(brewer.pal(3, "Purples")[2], brewer.pal(3, "Blues")[1:2])

## names and labels
names_models <- c("MV_PCA", "fPCA_off", "fPCA_gcv")
lables_models <- c("MV-PCA", "fPCA (no calibration)", "fPCA (gcv calibration)")

## resolution of the high resolution grid
n_nodes_HR_grid <- 1000


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Test: fPCA ----
cat.script_title(paste("Test:", TEST_SUITE))


## options ----
cat.section_title("Options")

## check arguments passed by terminal, set default if not present
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  RSTUDIO <- FALSE
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
domain_boundary <- generated_domain$domain_boundary
loadings_true_generator <- generated_domain$loadings_true_generator
mesh <- generated_domain$mesh

## locations
locations <- generate_locations(generated_domain, locs_eq_nodes)

## grid
grid <- spsample(domain_boundary, n_nodes_HR_grid, "regular")@coords

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

## load saved results
if (RUN$analysis) {
  source(paste("tests/", test_suite, "/templates/load_results.R", sep = ""))
}


### quantitative analysis ----
cat.subsection_title("Quantitative analysis")

## open a pdf where to save plots (quantitative analysis)
if (!RSTUDIO) {
  pdf(file = paste(path_images, name_test, "_quantitative.pdf", sep = ""))
}

## plots
if (RUN$quantitative_analysis) {
  source(paste("tests/", test_suite, "/templates/plot_quantitative_results.R", sep = ""))
}

## close pdf (quantitative analysis)
dev.off()


### qualitative analysis ----
cat.subsection_title("Qualitative analysis")

## open a pdf where to save plots (qualitative analysis)
if (!RSTUDIO) {
  pdf(file = paste(path_images, name_test, "_qualitative.pdf", sep = ""), width = 10, height = 7)
}

## plots
if (RUN$qualitative_analysis) {
  source(paste("tests/", test_suite, "/templates/plot_qualitative_results.R", sep = ""))
}

## close pdf (qualitative analysis)
dev.off()


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## remove options file
file.remove(paste(path_queue, file_options, sep = ""))

cat("\n\n")

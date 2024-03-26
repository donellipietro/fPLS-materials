# % %%%%%%%%%%%%%%%%%%%%%% %
# % % Test: Test Suite 1 % %
# % %%%%%%%%%%%%%%%%%%%%%% %

rm(list = ls())
graphics.off()


## global variables ----

test_suite <- "test_suite_1"
TEST_SUITE <- "Test Suite 1"


## prerequisite ----
source(paste("tests/", test_suite, "/utils/generate_options.R", sep = ""))


## libraries ----

## json
suppressMessages(library(jsonlite))

## data visualization
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gridExtra))


## sources ----
source("src/utils/directories.R")
source("src/utils/results_management.R")
source("src/utils/plots.R")


## paths ----
path_options <- paste("queue/", sep = "")
path_results <- paste("results/", test_suite, sep = "")
path_images <- paste("images/", test_suite, sep = "")
file_log <- "log.txt"


## options ----

## colors
colors <- c(brewer.pal(3, "Greens")[2], brewer.pal(3, "Blues")[2], brewer.pal(3, "Purples")[2])

## names and labels
names_models <- c("model1", "model2", "model2")
lables_models <- c("Model 1", "Model 2", "Model 3")


## load data ----

## check arguments passed by terminal
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  args[1] <- "test1"
}

## main test name
name_main_test <- args[1]
cat(paste("\nTest selected:", name_main_test, "\n"))
path_results <- paste(path_results, name_main_test, "/", sep = "")
path_images <- paste(path_images, name_main_test, "/", sep = "")

## generate options
generate_options(name_main_test, path_options)

## list of available tests
file_test_vect <- sort(list.files(path_options))

## open a pdf where to save the plots
pdf(paste(path_images, "time_complexity.pdf", sep = ""))

## room for solutions
columns_names <- c("Group", "n_locs", names_models)
empty_df <- data.frame(matrix(NaN, nrow = 0, ncol = 5))
colnames(empty_df) <- columns_names
times <- empty_df

for (file_test in file_test_vect) {
  ## load specs
  file_json <- paste(path_options, file_test, sep = "")
  parsed_json <- fromJSON(file_json)
  name_test <- parsed_json$test$name_test
  n_locs <- parsed_json$dimensions$n_locs
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
      list(
        n_locs = n_locs,
        fPCA_off = format_time(results_evaluation$fPCA_off$execution_time),
        fPCA_gcv = format_time(results_evaluation$fPCA_gcv$execution_time),
        SpatialPCA = format_time(results_evaluation$SpatialPCA$execution_time)
      ),
      columns_names
    )

    cat(paste("- Batch", i, "loaded\n"))
  }

  ## remove option file
  file.remove(file_json)
}

## plots
times <- times[, c("n_locs", "fPCA_off", "fPCA_gcv", "SpatialPCA")]
colnames(times) <- c("Group", "fPCA_off", "fPCA_gcv", "SpatialPCA")
boxplot <- plot.grouped_boxplots(
  times,
  values_name = "Time [seconds]",
  group_name = "S",
  subgroup_name = "Approches",
  subgroup_labels = c("fPCA (no calibration)", "fPCA (gcv)", "SpatialPCA"),
  subgroup_colors = colors,
  LEGEND = FALSE
) + standard_plot_settings() + ggtitle("Execution times")
plot <- plot.multiple_lines(
  times %>% group_by(Group) %>% summarize(fPCA_off = mean(fPCA_off), fPCA_gcv = mean(fPCA_gcv), SpatialPCA = mean(SpatialPCA)),
  values_name = "Time [seconds]",
  x_name = "S",
  subgroup_name = "Approches",
  subgroup_labels = c("fPCA (no calibration)", "fPCA (gcv)", "SpatialPCA"),
  subgroup_colors = colors,
  LOGLOG = TRUE,
  LEGEND = FALSE
) + standard_plot_settings() + ggtitle("Execution times - loglog scale")
print(boxplot)
print(plot)


dev.off()

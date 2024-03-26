# % %%%%%%%%%%%%%%%%%%%%%% %
# % % Test: Test Suite 1 % %
# % %%%%%%%%%%%%%%%%%%%%%% %

rm(list = ls())

## sources ----
source("src/utils/directories.R")
source(paste("tests/test_suite_1/utils/generate_options.R", sep = ""))


## libraries ----

## json
suppressMessages(library(jsonlite))


## directories ----
path_options <- paste("queue/", sep = "")
mkdir(c(path_options))


## options ----

## check arguments passed by terminal
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  args[1] <- "test1"
}

name_test_main <- args[1]

## main test name
name_main_test <- args[1]
cat(paste("\nTest selected:", name_test_main, "\n"))


## generation ----
generate_options(name_test_main, path_options)

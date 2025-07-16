# Network Imputation + Inference 10000

library(argparse)
library(minet)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

# Functions
random_replace <- function(dt, max_val) {
  dt[, (names(dt)) := lapply(.SD, function(x) {
    x[x == 0] <- runif(sum(x == 0), min = 0, max = max_val)
    return(x)
  })]
}

# Argument parser
parser <- ArgumentParser()
parser$add_argument("--input", required = TRUE, help = "Input directory")
parser$add_argument("--output", required = TRUE, help = "Output directory")
args <- parser$parse_args()

# Get list of files
file.list <- list.files(args$input, full.names = TRUE)
keep.list <- grep("10000", file.list, value = TRUE)
keep.list <- grep(".csv", keep.list, value = TRUE)

# Set up parallel backend
num_cores <- length(keep.list)  # Cores based on number of files in directory
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Run in parallel using foreach
#results <-
foreach(i = keep.list, .packages = c('data.table', 'minet', 'dplyr')) %dopar% {
  # Extract subject name
  subject_name <- tools::file_path_sans_ext(basename(i))
  cat(subject_name, "loading...\n")

  # Read and preprocess data
  patient <- read.csv(i)
  df <- patient %>% replace(is.na(.), 0)
  df <- as.matrix(df)

  # Imputation
  smallest_value <- (sort(unique(as.vector(df)))[2]) / 1e7
  df <- as.data.table(df)
  set.seed(7)
  imputed <- random_replace(df, smallest_value)
  imputed <- as.matrix(imputed)
  write.csv(imputed, file.path(args$output, paste0(subject_name, "_imputed.csv")), row.names = FALSE)

  # Build mutual information matrix
  mim_matrix <- build.mim(imputed, estimator = "spearman", disc = "equalfreq")
  write.csv(mim_matrix, file.path(args$output, paste0(subject_name, "_matrix.csv")), row.names = FALSE)

  # Infer networks
  network_1 <- aracne(mim_matrix)
  write.csv(network_1, file.path(args$output, paste0(subject_name, "_ARACNE.csv")), row.names = FALSE)

  network_2 <- clr(mim_matrix)
  write.csv(network_2, file.path(args$output, paste0(subject_name, "_CLR.csv")), row.names = FALSE)

  network_3 <- mrnet(mim_matrix)
  write.csv(network_3, file.path(args$output, paste0(subject_name, "_MRNET.csv")), row.names = FALSE)

  cat(subject_name, "complete!\n")
  return(i)
}

# Stop the cluster after the loop
stopCluster(cl)

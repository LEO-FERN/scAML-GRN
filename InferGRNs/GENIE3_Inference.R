# Network Inference GENIE3

library(argparse)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)
library(doRNG)
library(GENIE3)

# Argument parser
parser <- ArgumentParser()
parser$add_argument("--input", required = TRUE, help = "Input directory")
parser$add_argument("--output", required = TRUE, help = "Output directory")
args <- parser$parse_args()

# Get list of files
file.list <- list.files(args$input, full.names = TRUE)
keep.list <- grep("_imputed.csv", file.list, value = TRUE) # Select only imputed files

# Set up parallel backend
num_cores <- length(keep.list)*2
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Run in parallel using foreach
foreach(i = keep.list, .packages = c('data.table', 'dplyr', 'GENIE3', 'doRNG')) %dorng% {

  # Extract subject name from file name
  class.name <- gsub("10000_Genes_(.*)_imputed\\.csv$", "\\1", basename(i))
  cat(class.name, "loading...\n")

  # Read and preprocess data
  patient <- read.csv(i)
  df <- patient %>% replace(is.na(.), 0)  # Replace NA with zeros
  df <- as.matrix(df)

  # Add sample column as row names
  rownames(df) <- paste0("Sample_", seq_len(nrow(df)))

  # Transpose df
  df <- t(df)

  # For reproducibility of results
  set.seed(123)

  # Run GENIE3
  cat(class.name, "running GENIE3...\n")
  weightMatrix <- GENIE3(df, nCores = 2)

  # Write results
  output_file <- file.path(args$output, paste0(class.name, "_GENIE3.csv"))
  write.csv(weightMatrix, output_file, row.names = TRUE)
  cat(class.name, "completed: Results written to", output_file, "\n")

  return(output_file)
}

# Stop the cluster after the loop
stopCluster(cl)


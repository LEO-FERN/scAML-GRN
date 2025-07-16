library(netZooR)
library(data.table)
library(readr)
library(parallel)

# Function to process a single dataset
process_dataset <- function(file_path, output_dir) {
    # Read the CSV file
    data <- read_csv(file_path, col_names = TRUE)

    # Add sample index
    rownames(data) <- paste0("Sample", 1:nrow(data))

    # Transpose the dataframe
    transposed_data <- as.data.frame(t(data))

    # Ensure gene names remain as row names
    colnames(transposed_data) <- rownames(data)
    rownames(transposed_data) <- colnames(data)

    # Run LIONESS
    lioness_result <- lioness(transposed_data, network.inference.method = 'pearson', progress=TRUE)

    # Create output directory if it doesn't exist
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # Save each network as a CSV file
    for (i in 1:length(lioness_result)) {
        network <- lioness_result[[i]]
        output_file <- file.path(output_dir, paste0("network_", i, ".csv"))
        fwrite(as.data.table(network), output_file)
    }
}

# Function to process all files in a directory
process_directory <- function(input_dir, output_base_dir) {
    files <- list.files(input_dir, pattern = "_filtered.csv$", full.names = TRUE)
    for (file in files) {
        patient_id <- gsub("_filtered.csv", "", basename(file))
        output_dir <- file.path(output_base_dir, patient_id)
        process_dataset(file, output_dir)
        cat("Processed:", file, "\n")
    }
}

# Main execution
cell_types <- c("Dendritic", "Monocyte", "Progenitor")
base_input_dir <- "~/SingleCellData"
base_output_dir <- "~/SingleCellData/LIONESS_Output"

# Set up parallel processing
num_cores <- 3
cl <- makeCluster(num_cores)

# Export necessary functions and libraries to the cluster
clusterExport(cl, c("process_dataset", "process_directory", "base_input_dir", "base_output_dir"))
clusterEvalQ(cl, {
    library(netZooR)
    library(data.table)
    library(readr)
})

# Process directories in parallel
parLapply(cl, cell_types, function(cell_type) {
    input_dir <- file.path(base_input_dir, cell_type)
    output_dir <- file.path(base_output_dir, cell_type)
    cat("Processing", cell_type, "cells...\n")
    process_directory(input_dir, output_dir)
    cat("Finished processing", cell_type, "cells.\n\n")
})

# Stop the cluster
stopCluster(cl)

library(data.table)
library(minet)
library(GENIE3)
library(ggplot2)
library(parallel)

# Permutation Function
permute_dataset <- function(data) {
  permuted <- data
  for (col in colnames(data)) {
    permuted[[col]] <- sample(data[[col]])
  }
  return(permuted)
}

# Function to perform inference for all methods on a single permutation
infer_networks_for_permutation <- function(permuted_data, permutation_id, subject_name, output_dir) {
  results <- list()

  # Compute Mutual Information Matrix
  mim_matrix <- build.mim(as.matrix(permuted_data), estimator = "spearman", disc = "equalfreq")

  # ARACNE
  aracne_network <- aracne(mim_matrix)
  results$ARACNE <- as.vector(aracne_network[upper.tri(aracne_network, diag = FALSE)])

  # CLR
  clr_network <- clr(mim_matrix)
  results$CLR <- as.vector(clr_network[upper.tri(clr_network, diag = FALSE)])

  # MRNET
  mrnet_network <- mrnet(mim_matrix)
  results$MRNET <- as.vector(mrnet_network[upper.tri(mrnet_network, diag = FALSE)])

  # GENIE3 (transposed data for compatibility)
  transposed_data <- t(as.matrix(permuted_data))
  genie3_network <- GENIE3(transposed_data, nCores = 1)
  results$GENIE3 <- as.vector(genie3_network)

  # Save null distributions for each method
  for (method in names(results)) {
    null_file <- file.path(output_dir, paste0(subject_name, "_", method, "_perm_", permutation_id, ".csv"))
    write.csv(data.frame(edge_weights = results[[method]]), null_file, row.names = FALSE)
  }

  return(results)
}

# Main Function to Generate Null Distributions
generate_null_distributions <- function(file_path, output_dir, num_permutations = 100, num_cores = 25) {
  subject_name <- tools::file_path_sans_ext(basename(file_path))
  cat("Processing file:", file_path, "\n")

  # Load and preprocess data
  data <- fread(file_path)

  # Set up parallel cluster
  cluster <- makeCluster(num_cores)
  clusterExport(cluster, c("permute_dataset", "infer_networks_for_permutation", "build.mim",
                           "aracne", "clr", "mrnet", "GENIE3", "output_dir", "subject_name"),
                envir = environment())

  # Run permutations in parallel
  results_list <- parLapply(cluster, 1:num_permutations, function(perm_id) {
    cat("Processing Permutation:", perm_id, "\n")
    permuted_data <- permute_dataset(data)
    infer_networks_for_permutation(permuted_data, perm_id, subject_name, output_dir)
  })

  stopCluster(cluster)

  # Combine null distributions for each method
  combined_null_distributions <- list()
  for (method in c("ARACNE", "CLR", "MRNET", "GENIE3")) {
    combined_null_distributions[[method]] <- unlist(lapply(results_list, function(res) res[[method]]))

    # Save combined null distribution
    null_file <- file.path(output_dir, paste0(subject_name, "_", method, "_combined_null.csv"))
    write.csv(data.frame(edge_weights = combined_null_distributions[[method]]), null_file, row.names = FALSE)
  }

  return(combined_null_distributions)
}

# Visualize and Save Results
visualize_and_save_results <- function(null_distributions, output_dir, subject_name) {
  thresholds <- list()

  for (method in names(null_distributions)) {
    edge_weights <- null_distributions[[method]]

    # Plot null distribution
    p <- ggplot(data.frame(edge_weights = edge_weights), aes(x = edge_weights)) +
      geom_density(fill = "blue", alpha = 0.5) +
      labs(title = paste("Null Edge Weight Distribution -", method, "(", subject_name, ")"),
           x = "Edge Weight", y = "Density") +
      theme_minimal()

    # Save plot
    plot_file <- file.path(output_dir, paste0(subject_name, "_", method, "_null_distribution_plot.png"))
    ggsave(plot_file, p)

    # Calculate Threshold (e.g., 95th percentile)
    threshold <- quantile(edge_weights, 0.95)
    thresholds[[method]] <- threshold
  }

  # Save all thresholds to a single file
  thresholds_df <- data.frame(Method = names(thresholds), Threshold = unlist(thresholds))
  thresholds_file <- file.path(output_dir, paste0(subject_name, "_thresholds.csv"))
  write.csv(thresholds_df, thresholds_file, row.names = FALSE)

  return(thresholds)
}

# Main Script Execution
file_path <- "~/Data/Final_Dendritic_Net/AML556_imputed.csv"  # Provide the file path here
output_dir <- "~/Data/Final_Dendritic_Net"

# Parameters
num_permutations <- 50
num_cores <- 25

# Generate null distributions
null_distributions <- generate_null_distributions(file_path, output_dir, num_permutations, num_cores)

# Visualize and save thresholds
subject_name <- tools::file_path_sans_ext(basename(file_path))
thresholds <- visualize_and_save_results(null_distributions, output_dir, subject_name)

# Print thresholds
cat("Thresholds saved to:", file.path(output_dir, paste0(subject_name, "_thresholds.csv")), "\n")
print(thresholds)

library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(stringr)

# Define file paths
mapping_dir <- "~/EnrichmentData"
feature_importance_dir <- "~/ThesisCode/SingleCell/"
output_file <- "~/EnrichmentData/top_important_genes.csv"
output_dir <- "~/EnrichmentData"

# Function to process one cell type
process_cell_type <- function(cell_type) {
  # Read gene index mapping
  gene_mapping_file <- file.path(mapping_dir, paste0(cell_type, "_gene_index_mapping.csv"))
  gene_mapping <- read_csv(gene_mapping_file, show_col_types = FALSE)

  # Read feature importance
  feature_importance_file <- file.path(feature_importance_dir, paste0(cell_type, "_feature_importance.txt"))
  feature_importance <- read_table(feature_importance_file,
                                   skip = 2, # Skip the first two lines of metadata
                                   col_names = c("feature", "importance"),
                                   col_types = cols(feature = col_character(),
                                                    importance = col_double()))

  # Merge feature importance with gene mapping
  merged_data <- feature_importance %>%
    mutate(Index = as.numeric(str_extract(feature, "\\d+"))) %>%
    left_join(gene_mapping, by = "Index") %>%
    select(Gene1 = `Gene 1`, Gene2 = `Gene 2`, importance) %>%
    filter(!is.na(Gene1) & !is.na(Gene2))

  # Calculate number of features for top 0.05%
  #n_top_features <- nrow(merged_data) %/% 2000
  #print(n_top_features)

  # Get top 0.01% of features
  #top_features <- merged_data %>%
  #  top_n(n_top_features, importance)
  #print(top_features)

  # Get top 250 most important interactions
  top_features <- merged_data %>%
    slice_max(order_by = importance, n = 250)
  print(top_features)

  # Get unique genes from top features
  top_genes <- unique(c(top_features$Gene1, top_features$Gene2))
  top_genes <- top_genes[!is.na(top_genes)]
  #print(length(top_genes))

  # Select only the top 15 features
  top_15_features <- top_features %>% slice_max(order_by = importance, n = 15)

  # Save top 15 features to a separate file
  top_features_file <- file.path(output_dir, paste0(cell_type, "_top_15_features.csv"))
  write_csv(top_15_features, top_features_file)

  # Return results
  list(
    cell_type = cell_type,
    top_genes = paste(top_genes, collapse = ", ")
  )
}

# Process all cell types
cell_types <- c("Dendritic", "Monocyte", "Progenitor")
results <- lapply(cell_types, process_cell_type)

# Create a data frame with cell types and their important genes
important_genes_df <- map_df(results, ~data.frame(CellType = .x$cell_type,
                                                  ImportantGenes = .x$top_genes))

# Write the data frame to a CSV file in the EnrichmentData directory
write_csv(important_genes_df, output_file)

# Print a summary
cat("Top important genes have been saved to:", output_file, "\n")
for (result in results) {
  cat("\nCell Type:", result$cell_type, "\n")
  cat("Number of top genes:", length(str_split(result$top_genes, ", ")[[1]]), "\n")
}

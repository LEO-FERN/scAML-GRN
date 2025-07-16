# Load necessary libraries
library(data.table)
library(rlist)
library(dplyr)

# Function to merge datasets
merge.datasets <- function(files) {
  print("Read the datasets...")
  dataset.list <- list()
  class.labels <- c()

  for (x in files) {
    print(x)
    dataset <- fread(file = x)
    class.name <- gsub("^((AML\\d+[A|B]?)|(BM\\d+))_.*\\.csv$", "\\1", x)
    obs <- dim(dataset)[1]
    class.label <- rep(class.name, obs)
    class.labels <- append(class.labels, class.label)
    dataset.list <- rlist::list.append(dataset.list, dataset)
  }

  return(list(dataset.list, class.labels))
}

# Function to calculate the coefficient of variation
var.coeff <- function(x) {
  cv <- sd(x) / mean(x) * 100
  return(cv)
}

# Set the working directory
setwd('~/Data/Dendritic_Datasets/')

# Read files and merge datasets
file.list <- list.files()
keep.list <- grep("1000", file.list, value = TRUE, invert = TRUE) # invert = T, elements that do not match
keep.list <- grep(".csv", keep.list, value = TRUE)
merge.objects <- merge.datasets(files = keep.list)

# Expand the variables
sets.list <- merge.objects[[1]]
classes <- merge.objects[[2]]

# Merge all datasets into one
print("Merge the datasets...")
merged.dataset <- rbindlist(sets.list, fill = TRUE, use.names = TRUE)

# Preprocess the merged dataset
cols <- colnames(merged.dataset)
cols <- gsub('-', '.', cols)
colnames(merged.dataset) <- cols
merged.dataset[is.na(merged.dataset)] <- 0
merged.dataset <- as.matrix(merged.dataset[,-1])  # Remove metadata column

# Identify the top 10,000 most variable genes
filtering <- apply(merged.dataset, 2, var.coeff)
filter.sorted <- sort(filtering, decreasing = TRUE)
genes <- names(filter.sorted)[1:10000]

rm(merged.dataset, classes, cols, filtering)

# Process each patient dataset and write to CSV
for (file in keep.list) {
  # Read and preprocess the patient dataset
  AML <- as.matrix(read.csv(file))
  AML <- AML[, -1]  # Remove metadata column if present

  # Create a matrix with all 10,000 genes, fill missing genes with zero
  all_genes_matrix <- matrix(0, nrow = nrow(AML), ncol = length(genes))
  colnames(all_genes_matrix) <- genes
  rownames(all_genes_matrix) <- rownames(AML)

  # Populate the matrix with existing genes
  common_genes <- intersect(genes, colnames(AML))
  all_genes_matrix[, common_genes] <- AML[, common_genes, drop = FALSE]

  # Write the complete dataset to a new CSV file
  class.name <- gsub("^((AML\\d+[A|B]?)|(BM\\d+))_.*\\.csv$", "\\1", file)
  output_file <- paste0("10000_Genes_", class.name, ".csv")
  write.csv(all_genes_matrix, file = output_file, row.names = FALSE)
  print(output_file)
}

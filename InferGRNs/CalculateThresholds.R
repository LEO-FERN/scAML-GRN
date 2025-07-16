library(data.table)

setwd('~/Data/Final_Dendritic_Net/')

num_permutations <- 50
output_dir <- "."

combine_permutations <- function(technique, num_permutations) {
  combined_data <- vector("list", num_permutations)
  for (i in 1:num_permutations) {
    #filename <- sprintf("AML556_imputed_%s_perm_%d.csv", technique, i)
    filename <- sprintf("AML556_imputed_%s_perm_%d_SYM.csv", technique, i)
    data <- fread(filename)
    combined_data[[i]] <- data$edge_weights
  }
  return(unlist(combined_data))
}

process_technique <- function(technique, num_permutations) {
  cat("Processing", technique, "\n")

  edge_weights <- combine_permutations(technique, num_permutations)
  threshold <- quantile(edge_weights, 0.95)

  cat("Threshold for", technique, ":", threshold, "\n")

  return(data.frame(Technique = technique, Threshold = threshold))
}

main <- function() {
  #techniques <- c("ARACNE", "CLR", "GENIE3", "MRNET")
  techniques <- c("GENIE3")
  results <- data.frame(Technique = character(), Threshold = numeric())

  for (technique in techniques) {
    result <- process_technique(technique, num_permutations)
    results <- rbind(results, result)
  }

  write.csv(results, file.path(output_dir, "technique_thresholds_genie.csv"), row.names = FALSE)
  #cat("All thresholds have been calculated and saved to technique_thresholds.csv\n")
}

main()

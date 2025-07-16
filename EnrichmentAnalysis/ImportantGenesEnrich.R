library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(tidyr)

#install.packages("msigdbdf", repos = "https://igordot.r-universe.dev")

# Load your data
top_genes <- read.csv("~/EnrichmentData/top_important_genes.csv", stringsAsFactors = FALSE)

# Define gene sets
gene_sets <- list(
  "H"  = msigdbr(species = "Homo sapiens", collection = "H"),
  "C5" = msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "BP"),
  "C2" = msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "KEGG_MEDICUS"),
  "C6" = msigdbr(species = "Homo sapiens", collection = "C6")
)

# Choose your gene set (modify as needed)
selected_geneset_name <- "C6"  # Change to "H", "C5", etc.
selected_geneset <- gene_sets[[selected_geneset_name]]

# Prepare reference terms
reference.terms <- dplyr::select(selected_geneset, gs_name, gene_symbol)
reference.terms$gs_name <- gsub("KEGG_", "", reference.terms$gs_name)
reference.terms$gs_name <- gsub("_", " ", reference.terms$gs_name)

# Enrichment function
perform.enrichment <- function(gene.set, reference.terms) {
  gene.set <- unique(gene.set)
  enrichment <- enricher(
    gene = gene.set,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    TERM2GENE     = reference.terms
  )
  return(enrichment)
}

# Process each CellType for enrichment
results <- top_genes %>%
  mutate(ImportantGenes = strsplit(ImportantGenes, ",")) %>%
  group_by(CellType) %>%
  group_modify(~ {
    genes <- unlist(.x$ImportantGenes)
    genes <- trimws(genes)  # Remove whitespace

    # Perform enrichment analysis
    ego <- perform.enrichment(
      gene.set = genes,
      reference.terms = reference.terms
    )

    if (!is.null(ego) && nrow(ego@result) > 0) {
      sig_pathways <- ego@result %>%
        filter(p.adjust < 0.05) %>%
        pull(ID)
    } else {
      sig_pathways <- character(0)
    }

    tibble(Pathways = list(sig_pathways))
  }) %>%
  ungroup()

# Create unique pathway list
all_pathways <- unique(unlist(results$Pathways))
pathway_vec <- sort(all_pathways)

# Create binary matrix
binary_matrix <- results %>%
  mutate(PathwayVector = lapply(Pathways, function(x) {
    vec <- rep(0, length(pathway_vec))
    matches <- match(x, pathway_vec)
    matches <- matches[!is.na(matches)]
    vec[matches] <- 1
    vec
  })) %>%
  select(-Pathways) %>%
  unnest_wider(PathwayVector, names_sep = "") %>%
  rename_with(~ pathway_vec, starts_with("PathwayVector"))

# Final output
final_output <- binary_matrix %>%
  select(CellType, everything())

# Save results
write.csv(final_output, paste0("~/EnrichmentData/important_pathway_matrix_", selected_geneset_name, ".csv"), row.names = FALSE)

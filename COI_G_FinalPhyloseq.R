################################################################
# (COI-G) CREATE FINAL PHYLOSEQ OBJECT - DNA metabarcoding processing
################################################################
# Purpose: This script defines a function that takes the final LCA-assigned
#          taxonomy and the lulu-curated OTU table to create the final,
#          analysis-ready phyloseq object.
#
# Input:
#   - lca_output_path:       File path to the final taxonomy from galaxy-tool-lca.
#   - curated_phyloseq_path: File path to the lulu-curated phyloseq RDS object.
#   - final_phyloseq_path:   File path to save the final phyloseq RDS object.
#
# Output:
#   - The final, analysis-ready phyloseq object.
#   - The final phyloseq object is also saved to disk.
#
# Author: Paige Smallman, 2025
################################################################

create_final_phyloseq <- function(lca_output_path,
                                  curated_phyloseq_path,
                                  final_phyloseq_path) {

  # --- 1. Load Packages and Validate Inputs ---
  suppressPackageStartupMessages({
    require(phyloseq)
    require(dplyr)
    require(tidyr)
  })

  # Check if input files exist
  if (!file.exists(curated_phyloseq_path)) {
    stop("Curated phyloseq object not found at: ", curated_phyloseq_path)
  }
  if (!file.exists(lca_output_path) || file.info(lca_output_path)$size == 0) {
    stop("LCA output file not found or is empty at: ", lca_output_path,
         "\nPlease ensure the galaxy-tool-lca script was run successfully.")
  }

  # --- 2. Load Input Data ---
  message("\nLoading lulu-curated phyloseq object from: ", curated_phyloseq_path)
  ps_curated <- readRDS(curated_phyloseq_path)

  message("Loading LCA-assigned taxonomy from: ", lca_output_path)
  # The galaxy-lca tool output has headers that can be problematic for R.
  # We use `check.names = FALSE` to preserve them as is.
  tax_lca <- read.delim(lca_output_path, header = TRUE, sep = "\t", check.names = FALSE)

  # --- 3. Process Final Taxonomy Table ---
  message("Processing final taxonomy table...")

  # Keep only the first unique assignment for each OTU and set rownames
  tax_final_df <- tax_lca %>%
    distinct(`#Query`, .keep_all = TRUE) %>%
    column_to_rownames(var = "#Query")

  # Select and rename columns to standard taxonomic ranks
  tax_final_matrix <- tax_final_df %>%
    select(
      Kingdom = kingdom,
      Phylum = phylum,
      Class = class,
      Order = order,
      Family = family,
      Genus = genus,
      Species = species
    ) %>%
    as.matrix()

  # --- 4. Create Final Phyloseq Object ---
  message("Constructing the final phyloseq object...")
  
  # Get the list of OTUs that have a final taxonomic assignment
  final_otu_ids <- rownames(tax_final_matrix)
  
  # Prune the curated phyloseq object to keep only these OTUs.
  # This syncs the OTU table, sample data, and refseqs all at once.
  ps_final <- prune_taxa(final_otu_ids, ps_curated)
  
  # Now, replace the old (DADA2-based) taxonomy in the phyloseq object
  # with the new, final LCA-based taxonomy.
  # We subset the matrix by `taxa_names(ps_final)` to ensure the order is correct.
  tax_table(ps_final) <- tax_table(tax_final_matrix[taxa_names(ps_final), ])

  # --- 5. Save Final Object and Return ---
  message("\nFinal phyloseq object created with ", ntaxa(ps_final), " OTUs and ", nsamples(ps_final), " samples.")
  saveRDS(ps_final, file = final_phyloseq_path)
  message("Final phyloseq object saved to: ", final_phyloseq_path)

  print(ps_final)
  return(ps_final)
}
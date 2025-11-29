################################################################
# CREATE FINAL PHYLOSEQ OBJECT - DNA metabarcoding processing
################################################################

# This function takes the final LCA-assigned taxonomy and the lulu-curated
# OTU table to create the final, analysis-ready phyloseq object.

create_final_phyloseq <- function(lca_output_path,
                                  curated_phyloseq_path,
                                  final_phyloseq_path) {

  # Load required packages
  suppressPackageStartupMessages({
    require(phyloseq)
    require(dplyr)
    require(tidyr)
  })

  # --- 1. Load Input Data ---
  message("Loading lulu-curated phyloseq object from: ", curated_phyloseq_path)
  ps_curated <- readRDS(curated_phyloseq_path)

  message("Loading LCA-assigned taxonomy from: ", lca_output_path)
  # The galaxy-lca tool output has headers that get mangled by R's read.delim.
  # We read it in and clean up the column names.
  tax_lca <- read.delim(lca_output_path, header = TRUE, sep = "\t", check.names = FALSE)
  
  # --- 2. Process Taxonomy Table ---
  message("Processing final taxonomy table...")
  
  # Keep only the first unique assignment for each OTU
  tax_lca_unique <- tax_lca %>%
    distinct(`#Query`, .keep_all = TRUE)

  # Set OTU IDs as rownames
  rownames(tax_lca_unique) <- tax_lca_unique$`#Query`

  # Select and rename columns to standard taxonomic ranks
  tax_final <- tax_lca_unique %>%
    select(
      Kingdom = kingdom,
      Phylum = phylum,
      Class = class,
      Order = order,
      Family = family,
      Genus = genus,
      Species = species
    )

  # Convert to a matrix, which is required for the phyloseq tax_table
  tax_matrix <- as.matrix(tax_final)

  # --- 3. Synchronize OTU and Taxonomy Tables ---
  # Ensure that the OTUs in the OTU table and the new tax table match
  otu_table_final <- otu_table(ps_curated)
  
  # Find common OTUs between the final tax table and the OTU table
  common_otus <- intersect(rownames(otu_table_final), rownames(tax_matrix))
  
  message(paste("Found", length(common_otus), "common OTUs between OTU table and taxonomy."))

  # Prune both tables to include only the common OTUs
  otu_table_synced <- otu_table_final[common_otus, ]
  tax_matrix_synced <- tax_matrix[common_otus, ]

  # --- 4. Create Final Phyloseq Object ---
  message("Creating the final phyloseq object...")
  ps_final <- phyloseq(
    otu_table(otu_table_synced, taxa_are_rows = TRUE),
    tax_table(tax_matrix_synced),
    sample_data(ps_curated),
    refseq(ps_curated) # Keep the reference sequences
  )

  # --- 5. Save Final Object ---
  message("Saving final phyloseq object to: ", final_phyloseq_path)
  saveRDS(ps_final, file = final_phyloseq_path)

  print(ps_final)
  return(ps_final)
}
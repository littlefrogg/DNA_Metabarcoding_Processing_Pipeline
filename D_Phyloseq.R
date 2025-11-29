################################################################
# (D) CREATE INITIAL PHYLOSEQ OBJECT
################################################################
# Purpose: This script defines a function to combine the ASV table, taxonomy
#          table, and sample metadata into a phyloseq object.
#
# Input:
#   - otu_table:     The ASV table (matrix) from DADA2.
#   - tax_table:     The taxonomy table (matrix) from DADA2.
#   - metadata_path: File path to the sample metadata CSV file.
#   - output_path:   File path to save the final phyloseq RDS object.
#
# Output:
#   - A phyloseq object containing the combined data.
#   - The phyloseq object is also saved to disk at `output_path`.
#
# Author: Paige Smallman, 2025
################################################################

create_phyloseq <- function(otu_table,
                            tax_table,
                            metadata_path,
                            output_path) {

  # --- 1. Load Packages and Validate Inputs ---
  suppressPackageStartupMessages({
    require(phyloseq)
    require(tidyverse)
  })

  # Check that input objects are in the correct format
  if (!is.matrix(otu_table) || !is.matrix(tax_table)) {
    stop("Input 'otu_table' and 'tax_table' must be matrices.")
  }
  # Check that the metadata file exists
  if (!file.exists(metadata_path)) {
    stop("Metadata file not found at: ", metadata_path)
  }

  # --- 2. Load and Process Metadata ---
  message("\nLoading and processing metadata from: ", metadata_path)
  # Assuming the first column of the metadata CSV contains the sample IDs
  metadata <- read.csv(metadata_path, header = TRUE, row.names = 1)
  message("Metadata dimensions: ", paste(dim(metadata), collapse = " x "))

  # --- 3. Synchronize Samples between OTU Table and Metadata ---
  # Get sample IDs from both the OTU table and the metadata
  otu_samples <- colnames(otu_table)
  meta_samples <- rownames(metadata)

  # Find samples that are common to both
  common_samples <- intersect(otu_samples, meta_samples)
  
  # Report any mismatches
  missing_in_meta <- setdiff(otu_samples, meta_samples)
  missing_in_otu <- setdiff(meta_samples, otu_samples)

  if (length(missing_in_meta) > 0) {
    message("Warning: ", length(missing_in_meta), " samples from the OTU table are missing in the metadata and will be removed.")
  }
  if (length(missing_in_otu) > 0) {
    message("Warning: ", length(missing_in_otu), " samples from the metadata are missing in the OTU table and will be removed.")
  }

  # Subset both tables to include only the common samples
  otu_table_synced <- otu_table[, common_samples]
  metadata_synced <- metadata[common_samples, ]
  message("Retained ", length(common_samples), " matching samples for phyloseq object construction.")

  # --- 4. Construct and Save Phyloseq Object ---
  message("\nConstructing phyloseq object...")
  ps <- phyloseq(
    otu_table(otu_table_synced, taxa_are_rows = TRUE),
    sample_data(metadata_synced),
    tax_table(tax_table)
  )

  message("\nPhyloseq object created:")
  print(ps)

  # Save the final phyloseq object
  saveRDS(ps, file = output_path)
  message("Initial phyloseq object saved to: ", output_path)

  # Return the phyloseq object to the main script
  return(ps)
}

################################################################
# (C) INITIAL TAXONOMY ASSIGNMENT (DADA2)
################################################################
# Purpose: This script defines a function to assign taxonomy to ASVs using
#          the DADA2 `assignTaxonomy` function.
#
# Input:
#   - seqtab_nochim: The chimera-free ASV table (a matrix) from DADA2.
#   - tax_db_path:   File path to the DADA2-formatted taxonomy database.
#   - output_path:   File path to save the final taxonomy table RDS object.
#   - ncores:        Number of cores for parallel processing.
#   - seed:          A random seed for reproducibility.
#
# Output:
#   - A taxonomy table (matrix) with ASVs as rows and taxonomic ranks as columns.
#   - The taxonomy table is also saved to disk at `output_path`.
#
# Author: Paige Smallman, 2025
################################################################

assign_taxonomy <- function(seqtab_nochim,
                            tax_db_path,
                            output_path,
                            ncores = TRUE,
                            seed = 119) {

  # --- 1. Load Packages and Validate Inputs ---
  suppressPackageStartupMessages({
    require(dada2)
  })

  # Check that the taxonomy database file exists
  if (!file.exists(tax_db_path)) {
    stop("Taxonomy database not found at: ", tax_db_path)
  }
  # Check that the input is a matrix
  if (!is.matrix(seqtab_nochim)) {
    stop("Input 'seqtab_nochim' is not a matrix. Please provide a valid ASV table.")
  }

  # --- 2. Assign Taxonomy ---
  message("\nAssigning taxonomy using database: ", basename(tax_db_path))
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Time the execution for user feedback
  start_time <- Sys.time()
  
  tax_table <- assignTaxonomy(
    seqs = seqtab_nochim,
    refFasta = tax_db_path,
    multithread = ncores,
    tryRC = TRUE,  # Also check the reverse-complement for a match
    verbose = TRUE
  )
  
  end_time <- Sys.time()
  message("Taxonomy assignment completed in ", 
          round(difftime(end_time, start_time, units = "mins"), 1), 
          " minutes.")
  
  # --- 3. Save Output and Return ---
  # Save the final taxonomy table as an RDS object
  saveRDS(tax_table, file = output_path)
  message("Taxonomy table saved to: ", output_path)
  
  # Return the taxonomy table object to the main script
  return(tax_table)
}

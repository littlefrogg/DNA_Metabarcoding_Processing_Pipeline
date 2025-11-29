################################################################
# (COI-D) CURATE OTUS WITH LULU - DNA metabarcoding processing
################################################################
# Purpose: This script defines a function to curate an OTU table using the
#          LULU algorithm, which removes putatively erroneous OTUs based on
#          co-occurrence and similarity with more abundant "parent" OTUs.
#
# Input:
#   - ps_object:    A phyloseq object (must be OTU-clustered).
#   - matchlist_df: A data frame from VSEARCH containing the match list.
#   - output_path:  File path to save the final curated phyloseq RDS object.
#
# Output:
#   - A new, curated phyloseq object with erroneous OTUs removed.
#   - The curated phyloseq object is also saved to disk at `output_path`.
#
# Author: Paige Smallman, 2025
################################################################

curate_lulu <- function(ps_object, matchlist_df, output_path) {

  # --- 1. Load Packages and Validate Inputs ---
  suppressPackageStartupMessages({
    require(phyloseq)
    require(lulu)
    require(tidyverse)
  })

  # Check inputs
  if (!is(ps_object, "phyloseq")) {
    stop("Input 'ps_object' must be a phyloseq object.")
  }
  if (!is.data.frame(matchlist_df)) {
    stop("Input 'matchlist_df' must be a data frame.")
  }

  # --- 2. Prepare OTU Table and Run LULU ---
  # Extract the OTU table as a matrix, which is the required format for lulu
  otu_table_lulu <- as(otu_table(ps_object), "matrix")
  
  # Ensure taxa are rows for lulu processing
  if (taxa_are_rows(ps_object) == FALSE) {
    otu_table_lulu <- t(otu_table_lulu)
  }

  # Run the LULU curation algorithm
  message("\nStarting LULU curation...")
  curated_result <- lulu(otu_table_lulu, matchlist_df, minimum_match = 84)

  # --- 3. Filter Phyloseq Object ---
  # Get the list of OTU IDs that were kept after curation
  retained_otus <- rownames(curated_result$curated_table)
  
  message("LULU curation complete. ",
          length(retained_otus), " out of ",
          ntaxa(ps_object), " OTUs were retained.")
  message("A total of ", curated_result$discarded_count, " OTUs were discarded.")
  
  # Prune the original phyloseq object to keep only the curated OTUs
  ps_curated <- prune_taxa(retained_otus, ps_object)
  
  # IMPORTANT: Replace the OTU table in the pruned object with the lulu-curated table.
  # This ensures that read counts from discarded OTUs are correctly merged into their parent OTUs.
  otu_table(ps_curated) <- otu_table(curated_result$curated_table, taxa_are_rows = TRUE)

  # --- 4. Save Output and Return ---
  # Save the final curated phyloseq object
  saveRDS(ps_curated, file = output_path)
  message("Lulu-curated phyloseq object saved to: ", output_path)
  
  print(ps_curated)
  
  # Return the curated phyloseq object to the main script
  return(ps_curated)
}
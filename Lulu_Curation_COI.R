################################################################
# LULU OTU CURATION - DNA metabarcoding processing
################################################################

# This function curates an OTU table using the LULU algorithm, which removes
# putatively erroneous OTUs based on co-occurrence with more abundant OTUs.

curate_lulu <- function(ps_object, matchlist_df) {

  # Load required packages
  suppressPackageStartupMessages({
    require(phyloseq)
    require(lulu)
    require(tidyverse)
  })

  # Check inputs
  if (!is(ps_object, "phyloseq")) {
    stop("Input ps_object must be a phyloseq object.")
  }
  if (!is.data.frame(matchlist_df)) {
    stop("matchlist_df must be a data frame.")
  }

  # Prepare OTU table for LULU
  otu_table_lulu <- otu_table(ps_object) %>%
    as.data.frame()
    
  # DEFENSIVE CODING: Ensure only numeric columns are included
  # This step filters out any non-numeric columns that may have been added
  # in previous steps, such as OTU IDs or sequences.
  otu_table_lulu <- otu_table_lulu[, sapply(otu_table_lulu, is.numeric)]
  
  # Ensure row names are character for consistency
  rownames(otu_table_lulu) <- as.character(rownames(otu_table_lulu))

  # Run the LULU curation algorithm
  message("\n[", Sys.time(), "] Starting LULU curation...")
  curated_result <- lulu(otu_table_lulu, matchlist_df)

  # Extract the curated OTU table and a list of retained OTU IDs
  curated_table <- curated_result$curated_table
  retained_otus <- rownames(curated_table)

  message("LULU curation complete. ",
          nrow(curated_table), " out of ",
          nrow(otu_table_lulu), " OTUs were retained.")
  message("A total of ", curated_result$discarded_count, " OTUs were discarded.")
  
  # Filter the phyloseq object to keep only the curated OTUs
  curated_ps_object <- prune_taxa(retained_otus, ps_object)
  
  return(curated_ps_object)
}
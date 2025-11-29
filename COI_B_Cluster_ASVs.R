################################################################
# (COI-B) CLUSTER ASVs (DECIPHER) - DNA metabarcoding processing
################################################################
# Purpose: This script defines a function to cluster similar ASVs into OTUs
#          (at a 97% similarity threshold) using the DECIPHER package. It
#          then creates a new phyloseq object where the ASV table is collapsed
#          into an OTU table, and a consensus taxonomy and representative
#          sequence are generated for each OTU.
#
# Input:
#   - ps_object:   A phyloseq object (typically post-decontamination).
#   - output_path: File path to save the final clustered phyloseq RDS object.
#   - cutoff:      The similarity cutoff for clustering (e.g., 0.03 for 97%).
#   - processors:  Number of cores for parallel processing.
#
# Output:
#   - A new phyloseq object with OTUs instead of ASVs.
#   - The clustered phyloseq object is also saved to disk at `output_path`.
#
# Author: Paige Smallman, 2025
################################################################

cluster_coi_asvs <- function(ps_object,
                             output_path,
                             cutoff = 0.03,
                             processors = NULL) {

  # --- 1. Load Packages and Validate Inputs ---
  suppressPackageStartupMessages({
    require(phyloseq)
    require(DECIPHER)
    require(tidyverse)
  })

  # Extract sequences and validate
  dna <- refseq(ps_object)
  if (is.null(dna)) {
    stop("Phyloseq object is missing reference sequences (refseq), which are required for clustering.")
  }

  # --- 2. Cluster ASVs with DECIPHER ---
  message("\nStarting ASV clustering for ", ntaxa(ps_object), " ASVs...")
  clusters <- Clusterize(
    dna,
    cutoff = cutoff,
    processors = processors,
    verbose = TRUE
  )
  message("Clustering complete. Resulted in ", max(clusters$cluster), " OTUs.")

  # Create a mapping from ASV to OTU
  cluster_df <- tibble(
    ASV = names(dna),
    OTU = paste0("OTU_", clusters$cluster)
  )

  # --- 3. Create New OTU Table ---
  message("Creating new OTU table by merging ASV counts...")
  otu <- as(otu_table(ps_object), "matrix")
  if (taxa_are_rows(ps_object)) {
    otu <- t(otu)
  }
  # Ensure ASV order matches for safe aggregation
  cluster_df_sorted <- cluster_df[match(colnames(otu), cluster_df$ASV), ]
  otu_mat <- rowsum(t(otu), group = cluster_df_sorted$OTU)

  # --- 4. Generate Consensus Taxonomy for each OTU ---
  message("Generating consensus taxonomy for each OTU...")
  tax_df <- tax_table(ps_object) %>%
    as.data.frame() %>%
    rownames_to_column("ASV") %>%
    left_join(cluster_df, by = "ASV") %>%
    group_by(OTU) %>%
    # For each rank, find the most common non-NA assignment among the ASVs in the OTU
    summarize(across(everything(), ~{
      vals <- na.omit(.)
      if (length(vals) == 0) NA_character_ else names(which.max(table(vals)))
    })) %>%
    select(-ASV) %>% # Remove the now-irrelevant ASV column
    column_to_rownames("OTU") %>%
    as.matrix()

  # --- 5. Generate Representative Consensus Sequence for each OTU ---
  message("Generating representative sequence for each OTU...")
  # Create a DNAStringSet list where each element is the set of ASVs for one OTU
  dna_by_otu <- split(dna, cluster_df$OTU)
  # Generate a consensus sequence for each OTU
  rep_seqs <- DNAStringSet(sapply(dna_by_otu, function(x) as.character(ConsensusSequence(x)[[1]])))
  names(rep_seqs) <- names(dna_by_otu)
  # Ensure the order matches the new OTU table
  rep_seqs <- rep_seqs[rownames(otu_mat)]

  # --- 6. Build and Save Final OTU-based Phyloseq Object ---
  message("Constructing final OTU-based phyloseq object...")
  ps_otu <- phyloseq(
    otu_table(otu_mat, taxa_are_rows = TRUE),
    tax_table(tax_df),
    sample_data(ps_object),
    rep_seqs # Add the new representative sequences
  )

  # Save the final object
  saveRDS(ps_otu, output_path)
  message("\nSuccessfully saved clustered phyloseq object to: ", output_path)
  print(ps_otu)

  return(ps_otu)
}


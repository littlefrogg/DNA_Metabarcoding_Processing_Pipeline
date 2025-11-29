################################################################
# (COI-C) PREPARE FOR VSEARCH - DNA metabarcoding processing
################################################################
# Purpose: This script defines a function to extract the representative
#          sequences from an OTU-clustered phyloseq object and save them
#          as a FASTA file. This file is required as input for the VSEARCH
#          all-vs-all search used by LULU.
#
# Input:
#   - ps_object:         A phyloseq object (must be OTU-clustered).
#   - output_fasta_path: The full file path to save the output FASTA file.
#
# Output:
#   - A FASTA file containing the representative sequences of the OTUs.
#
# Author: Paige Smallman, 2025
################################################################

prepare_for_vsearch <- function(ps_object, output_fasta_path) {

  # --- 1. Load Packages and Validate Inputs ---
  suppressPackageStartupMessages({
    require(phyloseq)
    require(Biostrings)
  })

  # Extract reference sequences and check if they exist
  ref_seqs <- refseq(ps_object, errorIfNULL = FALSE)
  if (is.null(ref_seqs)) {
    stop("The input phyloseq object is missing reference sequences (refseq), which are required.")
  }

  # --- 2. Write Sequences to FASTA File ---
  message("\nPreparing FASTA file for VSEARCH...")
  
  # Use the robust writeXStringSet function to create the FASTA file
  Biostrings::writeXStringSet(ref_seqs, filepath = output_fasta_path)
  
  message("FASTA file of OTU representative sequences saved to: ", output_fasta_path)

  # --- 3. Return ---
  # This function's main purpose is to write a file, so it returns nothing.
  return(invisible(NULL))
}
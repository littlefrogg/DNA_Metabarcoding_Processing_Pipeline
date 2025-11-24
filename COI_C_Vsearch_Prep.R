################################################################
# CURATE OTUs AND PREPARE FOR VSEARCH - DNA metabarcoding
################################################################

# This function prepares an OTU table and FASTA file for external curation
# with tools like LULU or VSEARCH. It takes a phyloseq object as input.

curate_otus <- function(ps_object, output_dir, project_id) {

  # Load required packages
  suppressPackageStartupMessages({
    require(phyloseq)
    require(Biostrings)
    require(tidyverse)
  })

  # --- A: Prepare OTU table as a data frame (unchanged) ---
  otu_mat <- otu_table(ps_object) %>%
    as.data.frame() %>%
    t()

  otu_df <- otu_mat %>%
    as.data.frame() %>%
    rownames_to_column("OTU")

  # --- B: Prepare FASTA file for external tools (e.g., vsearch) ---
  message("\n[", Sys.time(), "] Preparing FASTA file for VSEARCH match list...")
  
  # Extract reference sequences
  ref_seqs <- refseq(ps_object)
  if (is.null(ref_seqs)) {
    stop("Phyloseq object missing reference sequences (refseq).")
  }

  # Create a data frame with the correct OTU IDs from your phyloseq object
  otu_seqs_df <- tibble(
    # Use the OTU names directly from the refseq object
    OTU = names(ref_seqs),
    Seq = as.character(ref_seqs)
  )

  # Create a FASTA-formatted character vector
  fa <- character(2 * nrow(otu_seqs_df))
  fa[c(TRUE, FALSE)] <- paste0(">", otu_seqs_df$OTU)
  fa[c(FALSE, TRUE)] <- otu_seqs_df$Seq

  fasta_file <- file.path(output_dir, paste0(project_id, "_otus.fasta"))
  writeLines(fa, con = fasta_file)
  message("Saved FASTA file to: ", fasta_file)

  # --- C: Save OTU table as a CSV for external use (unchanged) ---
  otu_file <- file.path(output_dir, paste0(project_id, "_otu_table.csv"))
  write.csv(otu_mat, file = otu_file, row.names = TRUE)
  message("Saved OTU table to: ", otu_file)

  message("\n[", Sys.time(), "] OTU curation step complete. Run VSEARCH with the generated FASTA file.")
  invisible(NULL)
}
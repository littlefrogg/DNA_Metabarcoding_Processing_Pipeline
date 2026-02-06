################################################################
# (COI-A) DECONTAMINATION - DNA metabarcoding processing
################################################################
# Purpose: This script defines a function to remove contaminant sequences
#          from a COI dataset. It applies several filtering steps:
#          1. Pre-filtering based on COI characteristics (length, N-bases).
#          2. Statistical contaminant removal using the decontam package.
#          3. Heuristic filtering for nuclear mitochondrial pseudogenes (NUMTs).
#          4. Filtering of low-read-count samples.
#
# Input:
#   - physeq:        A phyloseq object containing the initial ASV data.
#   - output_dir:    Path to the main output directory for the project.
#   - project_id:    A unique identifier for the project.
#   - control_col:   The column name in the metadata that identifies sample types.
#   - neg_controls:  A vector of names for negative control sample types.
#   - min_reads:     The minimum library size for a sample to be retained.
#   - max_n_ratio:   The maximum allowable ratio of 'N' bases in a sequence.
#
# Output:
#   - A list containing:
#     - phyloseq_clean: The final, decontaminated phyloseq object.
#     - contaminants: A data frame detailing which ASVs were flagged as contaminants.
#   - Diagnostic plots and a list of removed contaminants are saved to disk.
#
# Author: Paige Smallman, 2025
################################################################

run_COI_decontamination <- function(physeq,
                                    output_dir,
                                    project_id,
                                    control_col,
                                    neg_controls,
                                    min_reads,
                                    max_n_ratio,
                                    prevalence_threshold = 0.1,
                                    frame_shift_check = TRUE,
                                    manual_remove_taxa = NULL) {

  # --- 1. Load Packages and Validate Inputs ---
  suppressPackageStartupMessages({
    require(phyloseq)
    require(decontam)
    require(ggplot2)
    require(Biostrings)
    require(tidyverse)
  })

  # Validate that the input is a phyloseq object
  if (!is(physeq, "phyloseq")) {
    stop("Input 'physeq' must be a valid phyloseq object.")
  }
  # Check for reference sequences, which are required for COI-specific filters
  if (is.null(refseq(physeq, errorIfNULL = FALSE))) {
    stop("The phyloseq object must contain reference sequences (refseq) for this function.")
  }

  # --- 2. Extract Sequences and Apply COI Pre-filters ---
  message("\nApplying COI-specific pre-filters...")
  asv_seqs <- refseq(physeq)
  asv_seqs_all <- asv_seqs # Keep a copy of all original sequences for later
  
  # Debug: Print initial dimensions
  message("DEBUG: Initial phyloseq dimensions: ", paste(dim(otu_table(physeq)), collapse = " x "))

  # Filter 1: Remove sequences with a high proportion of 'N' bases
  n_content <- letterFrequency(asv_seqs, "N", as.prob = TRUE)
  keep_n <- n_content <= max_n_ratio

  # Filter 2: Check for frameshifts by ensuring sequence length is a multiple of 3
  if (frame_shift_check) {
    keep_frame <- width(asv_seqs) %% 3 == 0
  } else {
    keep_frame <- rep(TRUE, ntaxa(physeq))
  }

  # Combine filters and prune the phyloseq object
  keep_taxa_prefilter <- as.vector(keep_n & keep_frame)
  names(keep_taxa_prefilter) <- taxa_names(physeq)
  ps <- prune_taxa(keep_taxa_prefilter, physeq)
  asv_seqs <- refseq(ps) # Update asv_seqs to match the pruned object
  message(ntaxa(physeq) - ntaxa(ps), " ASVs removed by COI pre-filtering (N-content/frameshift).")
  
  # Debug: Print dimensions after pre-filtering
  message("DEBUG: Dimensions after pre-filtering: ", paste(dim(otu_table(ps)), collapse = " x "))

  # --- 3. Identify Contaminants with Decontam ---
  message("\nIdentifying contaminants with decontam (prevalence method)...")
  sample_data(ps)$is_neg <- sample_data(ps)[[control_col]] %in% neg_controls
  contam_df <- isContaminant(ps, method = "prevalence", neg = "is_neg", threshold = prevalence_threshold)

  # --- 4. Identify Potential NUMTs ---
  message("Applying NUMT detection heuristics (GC content and stop codons)...")
  # Heuristic 1: Identify outliers in GC content
  gc_content <- letterFrequency(asv_seqs, "GC", as.prob = TRUE)[, 1]
  gc_outliers <- gc_content < quantile(gc_content, 0.01) | gc_content > quantile(gc_content, 0.99)

  # Heuristic 2: Check for internal stop codons
  # Helper function to safely translate DNA, handling potential errors
  safe_translate <- function(dna) {
    sapply(dna, function(s) {
      tryCatch({
        # Use Invertebrate Mitochondrial Code (SGC5) to avoid false positives (e.g. TGA is Trp, not Stop)
        as.character(Biostrings::translate(s, genetic.code = Biostrings::getGeneticCode("SGC5")))
      }, error = function(e) "")
    })
  }
  aa_seqs <- safe_translate(asv_seqs)
  has_stop <- grepl("\\*", aa_seqs)
  
  # Ensure vectors are named to match phyloseq taxa
  names(gc_outliers) <- taxa_names(ps)
  names(has_stop) <- taxa_names(ps)

  # Combine NUMT indicators and add to the contaminant dataframe
  numt_candidates <- gc_outliers | has_stop
  
  # Initialize columns to avoid NA issues
  contam_df$numt_candidate <- FALSE
  if(all(rownames(contam_df) %in% names(numt_candidates))) {
      contam_df$numt_candidate <- numt_candidates[rownames(contam_df)]
  } else {
      # Fallback: assume order matches if names don't (shouldn't happen)
      contam_df$numt_candidate <- numt_candidates
  }
  
  # Handle any NAs that might have crept in
  contam_df$contaminant[is.na(contam_df$contaminant)] <- FALSE
  contam_df$numt_candidate[is.na(contam_df$numt_candidate)] <- FALSE
  
  # Flag any ASV that is either a decontam contaminant OR a NUMT candidate for removal
  contam_df$remove <- contam_df$contaminant | contam_df$numt_candidate

  # --- DEBUG: Print detailed removal stats ---
  message("DEBUG: Decontam identified ", sum(contam_df$contaminant, na.rm=TRUE), " contaminants.")
  message("DEBUG: NUMT check identified ", sum(contam_df$numt_candidate, na.rm=TRUE), " candidates.")
  message("DEBUG:   - GC outliers: ", sum(gc_outliers, na.rm=TRUE))
  message("DEBUG:   - Stop codons: ", sum(has_stop, na.rm=TRUE))
  message("DEBUG: Total ASVs to remove: ", sum(contam_df$remove, na.rm=TRUE))
  message("DEBUG: Total ASVs remaining: ", sum(!contam_df$remove, na.rm=TRUE))
  # -------------------------------------------

  # --- 5. Generate Diagnostic Plots ---
  contam_plot_data <- data.frame(p = contam_df$p, Contaminant = contam_df$contaminant)
  contam_hist <- ggplot(contam_plot_data, aes(x = p)) +
    geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") +
    labs(title = "Decontam P-value Distribution", x = "P-Value", y = "Frequency") +
    theme_minimal()
  ggsave(file.path(output_dir, "decontam_p_value_hist.png"), contam_hist, width = 8, height = 6)

  # --- 6. Filter Phyloseq Object ---
  message("\nRemoving contaminants and NUMTs...")
  
  # Keep only taxa that are NOT marked for removal
  keep_taxa_final <- rownames(contam_df)[!contam_df$remove]
  
  # Also remove manually specified taxa if any
  if (!is.null(manual_remove_taxa)) {
    keep_taxa_final <- setdiff(keep_taxa_final, manual_remove_taxa)
  }
  
  if (length(keep_taxa_final) == 0) {
      stop("DEBUG ERROR: All taxa were filtered out! Check the debug messages above to see which filter is too strict.")
  }
  
  physeq_clean <- prune_taxa(keep_taxa_final, ps)
  
  # Debug: Print dimensions after contaminant removal
  message("DEBUG: Dimensions after contaminant removal: ", paste(dim(otu_table(physeq_clean)), collapse = " x "))
  
  # Filter samples with low read counts
  physeq_clean <- prune_samples(sample_sums(physeq_clean) >= min_reads, physeq_clean)
  
  # Debug: Print dimensions after sample filtering
  message("DEBUG: Dimensions after sample filtering: ", paste(dim(otu_table(physeq_clean)), collapse = " x "))
  
  # Remove any taxa that are no longer present in any sample
  physeq_clean <- prune_taxa(taxa_sums(physeq_clean) > 0, physeq_clean)
  
  # Debug: Print final dimensions
  message("DEBUG: Final dimensions: ", paste(dim(otu_table(physeq_clean)), collapse = " x "))
  taxa_to_keep <- rownames(contam_df[contam_df$remove == FALSE, ])
  ps_clean <- prune_taxa(taxa_to_keep, ps)
  message(nrow(contam_df) - length(taxa_to_keep), " ASVs removed as contaminants or NUMTs.")

  # Remove samples with fewer than `min_reads`
  ps_clean <- prune_samples(sample_sums(ps_clean) >= min_reads, ps_clean)
  message(nsamples(ps) - nsamples(ps_clean), " samples removed due to low read counts (<", min_reads, ").")

  # Remove any ASVs that now have zero reads across all samples
  ps_clean <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)

  # Optional manual removal of specific taxa
  if (!is.null(manual_remove_taxa)) {
    taxa_to_remove_manual <- intersect(manual_remove_taxa, taxa_names(ps_clean))
    ps_clean <- prune_taxa(!taxa_names(ps_clean) %in% taxa_to_remove_manual, ps_clean)
    message("Manually removed ", length(taxa_to_remove_manual), " specified taxa.")
  }

  # --- 7. Save Outputs and Return ---
  # Save the list of removed contaminants for inspection
  removed_taxa_df <- contam_df[contam_df$remove == TRUE, ]
  removed_seqs <- asv_seqs_all[rownames(removed_taxa_df)]
  removed_data <- data.frame(
    ASV = rownames(removed_taxa_df),
    Sequence = as.character(removed_seqs),
    Decontam_Contaminant = removed_taxa_df$contaminant,
    NUMT_Candidate = removed_taxa_df$numt_candidate,
    Decontam_P_Value = removed_taxa_df$p
  )
  write.csv(removed_data, file.path(output_dir, "coi_removed_contaminants.csv"), row.names = FALSE)

  message("\nDecontamination complete. Final dimensions:")
  print(ps_clean)

  # Return the cleaned phyloseq object and the contaminant data frame
  return(list(
    phyloseq_clean = ps_clean,
    contaminants = contam_df
  ))
}

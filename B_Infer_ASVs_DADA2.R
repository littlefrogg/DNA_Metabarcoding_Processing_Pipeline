################################################################
# (B) INFER ASVs WITH DADA2
################################################################
# Purpose: This script defines a function to filter and trim reads, learn
#          error rates, infer ASVs, merge paired-end reads, and remove
#          chimeras using the DADA2 pipeline.
#
# Input:
#   - fnFs, fnRs: File paths for forward and reverse reads.
#   - sample.names: A character vector of sample names.
#   - pathinput: Path to the directory containing trimmed FASTQ files.
#   - pathoutput_tracktable: Path to save the read tracking table.
#   - pathfigures: Path to save error model plots.
#   - trunc_params: A list with DADA2 trimming and filtering parameters.
#   - ncores: Number of cores for parallel processing.
#   - use_nova_seq: Boolean flag to use NovaSeq-specific error model.
#   - use_coi_seqtab: Boolean flag to filter ASV table by COI length.
#
# Output:
#   - A list containing:
#     - seqtab_nochim: The final, chimera-free ASV table.
#     - track: A table tracking read counts through the pipeline.
#   - The read tracking table and error plots are saved to disk.
#
# Author: Paige Smallman, 2025
################################################################

run_dada2_processing <- function(fnFs,
                                 fnRs,
                                 sample.names,
                                 pathinput,
                                 pathoutput_tracktable,
                                 pathfigures,
                                 trunc_params,
                                 ncores = TRUE,
                                 use_nova_seq = FALSE,
                                 use_coi_seqtab = FALSE) {

  # --- 1. Setup and Parameter Validation ---
  # Load required packages
  suppressPackageStartupMessages({
    require(dada2)
    require(ggplot2)
  })

  # Create a subdirectory for filtered reads
  filtered_path <- file.path(pathinput, "filtered")
  dir.create(filtered_path, showWarnings = FALSE)

  # Validate that trunc_params are correctly formatted
  if (length(trunc_params$truncLen) != 2 || length(trunc_params$maxEE) != 2) {
    stop("Parameter error: 'truncLen' and 'maxEE' must both be vectors of length 2.")
  }
  message("\nStarting DADA2 processing with parameters:")
  message(paste(capture.output(print(trunc_params)), collapse = "\n"))

  # Verify files exist and are not empty
  # This filtering prevents parallel processing crashes on empty/tiny files (<25KB ~50 lines)
  min_file_size <- 25000
  size_check_passed <- file.size(fnFs) >= min_file_size & file.size(fnRs) >= min_file_size

  # Diagnostic: print and save raw file sizes for inspection
  raw_sizes <- data.frame(
    sample = sample.names,
    fnF = fnFs,
    fnR = fnRs,
    size_F = file.size(fnFs),
    size_R = file.size(fnRs),
    pass_size_filter = size_check_passed
  )
  message("\nRaw file size summary (first 10 shown):")
  print(head(raw_sizes, 10))
  write.csv(raw_sizes, file.path(pathinput, "filtered", "raw_file_size_summary.csv"), row.names = FALSE)

  if (any(!size_check_passed)) {
    dropped <- sample.names[!size_check_passed]
    message("\nWARNING: Excluding ", length(dropped), " samples with file size < ", min_file_size, " bytes (too few reads to process).")
    message("Excluded samples: ", paste(dropped, collapse = ", "))
    # Update vectors to keep only passed files
    fnFs <- fnFs[size_check_passed]
    fnRs <- fnRs[size_check_passed]
    sample.names <- sample.names[size_check_passed]
    if(length(fnFs) == 0) stop("No valid samples remained after filtering for minimum file size!")
  }

  # --- 2. Filter and Trim Reads ---
  # Define paths for filtered output files
  filtFs <- file.path(filtered_path, paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(filtered_path, paste0(sample.names, "_R_filt.fastq.gz"))
  names(filtFs) <- names(filtRs) <- sample.names

  # Perform filtering and trimming
  filter_stats <- filterAndTrim(
    fwd = fnFs, filt = filtFs,
    rev = fnRs, filt.rev = filtRs,
    truncLen = trunc_params$truncLen,
    maxN = 0,
    maxEE = trunc_params$maxEE,
    truncQ = trunc_params$truncQ,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = ncores
  )

  # Diagnostic: print and save filter stats for inspection
  filter_stats_df <- as.data.frame(filter_stats)
  filter_stats_df$sample <- sample.names
  filter_stats_df$fnF <- filtFs
  filter_stats_df$fnR <- filtRs
  message("\nFilter stats summary (first 10 shown):")
  print(head(filter_stats_df, 10))
  write.csv(filter_stats_df, file.path(filtered_path, "filter_stats_summary.csv"), row.names = FALSE)


  # Remove any samples that had zero reads pass filtering or resulted in near-empty files
  keep <- file.exists(filtFs) & file.exists(filtRs)
  if (any(!keep)) {
    message("\nWARNING: Excluding ", sum(!keep), " samples with missing filtered files.")
  }
  filtFs <- filtFs[keep]
  filtRs <- filtRs[keep]
  sample.names <- sample.names[keep]
  filter_stats <- filter_stats[keep, , drop=FALSE]
  
  if (length(filtFs) == 0) {
    stop("No samples remaining after filtering. Check input files and parameters.")
  }

  # Check for files that ended up with 0 reads after filtering (very small gzipped files)
  min_filtered_size <- 50000
  valid_indices <- file.size(filtFs) > min_filtered_size & file.size(filtRs) > min_filtered_size
  if (any(!valid_indices)) {
    dropped_final <- sample.names[!valid_indices]
    message("\nWARNING: Additional ", length(dropped_final), " samples dropped after quality filtering resulted in near-empty files.")
    message("Dropped post-filter: ", paste(dropped_final, collapse=", "))
    filtFs <- filtFs[valid_indices]
    filtRs <- filtRs[valid_indices]
    sample.names <- sample.names[valid_indices]
    filter_stats <- filter_stats[valid_indices, , drop=FALSE]
  }

  # Print vector lengths for debugging
  message(sprintf("\nSamples after filtering: %d", length(sample.names)))
  message(sprintf("filtFs: %d, filtRs: %d, filter_stats rows: %d", length(filtFs), length(filtRs), nrow(filter_stats)))
  if (!(length(filtFs) == length(filtRs) && length(filtFs) == length(sample.names) && length(filtFs) == nrow(filter_stats))) {
    stop("Vector length mismatch after filtering! Please check filtering logic.")
  }

  # --- 3. Learn Error Rates ---

  message("\nLearning error models...")
  if (use_nova_seq) {
    message("Using error model for binned NovaSeq quality scores.")
    errF <- learnErrors(filtFs, multithread = ncores, nbins = 16)
    errR <- learnErrors(filtRs, multithread = ncores, nbins = 16)
  } else {
    message("Using standard error model for Illumina data.")
    errF <- learnErrors(filtFs, multithread = ncores, randomize = TRUE)
    errR <- learnErrors(filtRs, multithread = ncores, randomize = TRUE)
  }

  # Save error model plots
  ggsave(file.path(pathfigures, "01_error_model_forward.png"), plotErrors(errF, nominalQ = TRUE))
  ggsave(file.path(pathfigures, "01_error_model_reverse.png"), plotErrors(errR, nominalQ = TRUE))

  # --- 4. Infer ASVs ---
  # Dereplicate reads
  derepFs <- derepFastq(filtFs, verbose = TRUE)
  derepRs <- derepFastq(filtRs, verbose = TRUE)
  names(derepFs) <- names(derepRs) <- names(filtFs)

  # Infer ASVs (the core DADA2 algorithm)
  dadaFs <- dada(derepFs, err = errF, multithread = ncores)
  dadaRs <- dada(derepRs, err = errR, multithread = ncores)

  # --- 5. Merge Reads and Construct ASV Table ---
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
  seqtab <- makeSequenceTable(mergers)
  message("\nInitial ASV table dimensions: ", paste(dim(seqtab), collapse = " x "))

  # Optional: Filter ASV table by expected length for COI
  if (use_coi_seqtab) {
    message("Filtering ASVs to keep only those with length 310–316 bp.")
    seqtab <- seqtab[, nchar(colnames(seqtab)) %in% 310:316]
    message("Dimensions after length filtering: ", paste(dim(seqtab), collapse = " x "))
  }

  # --- 6. Remove Chimeras ---
  seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = ncores, verbose = TRUE)
  message("Post-chimera removal dimensions: ", paste(dim(seqtab_nochim), collapse = " x "))

  # --- 7. Track Reads and Finalize Output ---
  getN <- function(x) sum(getUniques(x))
  
  # Start with the stats from the filtering step
  track <- as.data.frame(filter_stats)
  
  # Add columns one by one, matching by sample name
  track$denoisedF <- sapply(names(filtFs), function(sn) getN(dadaFs[[sn]]))
  track$denoisedR <- sapply(names(filtFs), function(sn) getN(dadaRs[[sn]]))
  track$merged <- sapply(names(filtFs), function(sn) getN(mergers[[sn]]))
  
  # For the non-chimera counts, handle samples that might have dropped out
  track$nonchim <- 0 # Default to 0
  present_samples <- rownames(track)[rownames(track) %in% rownames(seqtab_nochim)]
  track[present_samples, "nonchim"] <- rowSums(seqtab_nochim[present_samples, , drop = FALSE])

  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  
  # Save the tracking table
  write.table(track, pathoutput_tracktable, sep = "\t", quote = FALSE, col.names = NA)
  message("\nRead tracking table saved to: ", pathoutput_tracktable)

  # Return the key results to the main script
  return(list(
    seqtab_nochim = seqtab_nochim,
    track = track
  ))
}
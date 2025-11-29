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

  # Remove any samples that had zero reads pass filtering
  keep <- file.exists(filtFs) & file.exists(filtRs)
  filtFs <- filtFs[keep]
  filtRs <- filtRs[keep]
  if (length(filtFs) == 0) {
    stop("No samples remaining after filtering. Check input files and parameters.")
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
  track <- cbind(
    filter_stats[keep, ],
    sapply(dadaFs[names(filtFs)], getN),
    sapply(dadaRs[names(filtFs)], getN),
    sapply(mergers[names(filtFs)], getN),
    rowSums(seqtab_nochim[rownames(filter_stats[keep,]), ])
  )
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
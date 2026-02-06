################################################################
# (A) INITIAL QUALITY ASSESSMENT
################################################################
# Purpose: This script defines a function to find paired-end sequence reads,
#          ensure they are valid, and generate quality profile plots to help
#          determine trimming parameters for the DADA2 step.
#
# Input:
#   - pathinput:   Path to the directory containing trimmed FASTQ files.
#   - pathfigures: Path to the directory where output plots will be saved.
#
# Output:
#   - A list containing:
#     - fnFs: A character vector of forward read file paths.
#     - fnRs: A character vector of reverse read file paths.
#     - sample.names: A character vector of the sample names.
#   - PNG files of the quality plots are saved to `pathfigures`.
#
# Author: Paige Smallman, 2025
################################################################

generate_quality_plots <- function(pathinput, pathfigures) {
  
  # --- 1. Find Sequence Files with Automatic Pattern Detection ---
  # Define possible filename patterns for forward and reverse reads.
  patterns <- list(
    R1 = c("_R1_001.trimmed.fastq$", "_R1_trimmed.fastq$"),
    R2 = c("_R2_001.trimmed.fastq$", "_R2_trimmed.fastq$")
  )
  
  # Loop through patterns to find which one matches the files in the directory.
  fnFs <- NULL
  fnRs <- NULL
  for (i in seq_along(patterns$R1)) {
    fnFs_try <- sort(list.files(pathinput, pattern = patterns$R1[i], full.names = TRUE))
    fnRs_try <- sort(list.files(pathinput, pattern = patterns$R2[i], full.names = TRUE))
    
    # If files are found with this pattern, use them and stop searching.
    if (length(fnFs_try) > 0 && length(fnRs_try) > 0) {
      fnFs <- fnFs_try
      fnRs <- fnRs_try
      message("Found sequence files with pattern: ", patterns$R1[i], " and ", patterns$R2[i])
      break
    }
  }
  
  # If no matching files were found after checking all patterns, stop with an error.
  if (is.null(fnFs) || length(fnFs) == 0) {
    stop("No matching FASTQ files found in the input directory. Checked for patterns ending in '_R1_001.trimmed.fastq' and '_R1_trimmed.fastq'.")
  }
  
  # --- 2. Filter and Pair Files ---
  # Remove any files that are likely empty or have too few reads.
  min_file_size <- 1000 # 1KB
  
  # Check for small files and warn
  small_Fs <- fnFs[file.size(fnFs) < min_file_size]
  if(length(small_Fs) > 0) {
    message("Filtering out ", length(small_Fs), " forward files smaller than ", min_file_size, " bytes.")
    message("Dropped files: ", paste(basename(small_Fs), collapse=", "))
  }
  
  fnFs <- fnFs[file.size(fnFs) >= min_file_size]
  fnRs <- fnRs[file.size(fnRs) >= min_file_size]
  
  # Extract sample names by removing the suffix. This regex handles both patterns and separators.
  get_sample_name <- function(x) sub("_R[12](_001)?[._]trimmed\\.fastq$", "", basename(x))
  sample.names.F <- sapply(fnFs, get_sample_name)
  sample.names.R <- sapply(fnRs, get_sample_name)
  
  # Keep only the files that have a matching pair.
  common.samples <- intersect(sample.names.F, sample.names.R)
  fnFs <- fnFs[sample.names.F %in% common.samples]
  fnRs <- fnRs[sample.names.R %in% common.samples]
  sample.names <- common.samples
  
  # --- 3. Generate and Save Quality Plots ---
  # Validate files for plotting to prevent errors with corrupted/low-quality files
  validate_for_plotting <- function(files) {
    valid_files <- character(0)
    # Check max 20 smallest files + 20 random others to save time, or just all?
    # Checking all is safer.
    message("Validating valid FASTQ structure for ", length(files), " files (this may take a moment)...")
    
    for (f in files) {
      skip <- FALSE
      tryCatch({
        # Try plotting just 100 reads to check integrity
        invisible(dada2::plotQualityProfile(f, n = 100))
      }, error = function(e) {
        message("Skipping file for plotting (error): ", basename(f), " - ", e$message)
        skip <<- TRUE
      }, warning = function(w) {
         # Warnings like "no non-missing args" are what we want to skip
         if(grepl("no non-missing arguments", w$message) || grepl("minScore", w$message)) {
             message("Skipping file for plotting (warning): ", basename(f), " - ", w$message)
             skip <<- TRUE
         }
      })
      if (!skip) valid_files <- c(valid_files, f)
    }
    return(valid_files)
  }

  fnFs_plot <- validate_for_plotting(fnFs)
  fnRs_plot <- validate_for_plotting(fnRs)
  
  if(length(fnFs_plot) == 0) stop("No valid forward files found for plotting!")
  
  # Create the quality profile plots.
  qp_fwd <- plotQualityProfile(fnFs_plot, aggregate = TRUE) + 
    ggtitle("Forward Reads Quality Profile")
  qp_rev <- plotQualityProfile(fnRs_plot, aggregate = TRUE) + 
    ggtitle("Reverse Reads Quality Profile") 
  
  # Save the plots to the specified figures directory.
  dir.create(pathfigures, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(pathfigures, "00_quality_forward.png"), qp_fwd, width = 10, height = 6)
  ggsave(file.path(pathfigures, "00_quality_reverse.png"), qp_rev, width = 10, height = 6)
  
  # --- 4. Return File Paths and Sample Names ---
  # Return the validated file lists for the next step in the pipeline.
  return(list(fnFs = fnFs, fnRs = fnRs, sample.names = sample.names))
}

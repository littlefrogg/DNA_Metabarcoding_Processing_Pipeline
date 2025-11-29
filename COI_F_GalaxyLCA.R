################################################################
# (COI-F) PREPARE FOR GALAXY-LCA - DNA metabarcoding processing
################################################################
# Purpose: This script defines a function to format the BLAST output to be
#          compatible with the galaxy-tool-lca python script for Lowest
#          Common Ancestor (LCA) taxonomic assignment.
#
# Input:
#   - blast_output_path: File path to the tabular BLAST output file.
#   - lca_input_path:    File path to save the formatted output for the LCA tool.
#   - db_source_name:    A string identifying the source database (e.g., "MIDORI2").
#
# Output:
#   - A tab-delimited text file formatted for the galaxy-tool-lca script.
#   - The path to the created file is returned.
#
# Author: Paige Smallman, 2025
################################################################

format_blast_for_lca <- function(blast_output_path, lca_input_path, db_source_name = "MIDORI2") {

  # --- 1. Load Packages and Validate Inputs ---
  suppressPackageStartupMessages({
    require(dplyr)
    require(tidyr)
    require(stringr)
  })

  # Check if the BLAST output file exists and is not empty
  if (!file.exists(blast_output_path) || file.info(blast_output_path)$size == 0) {
    stop("BLAST output file not found or is empty at: ", blast_output_path,
         "\nPlease ensure the BLAST step completed successfully.")
  }

  # --- 2. Read and Format BLAST Output ---
  # Define column names based on the custom BLAST output format from the .sh script
  blast_cols <- c(
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs"
  )

  # Read BLAST output
  message("\nReading BLAST output from: ", blast_output_path)
  blast_out <- read.delim(blast_output_path, header = FALSE, col.names = blast_cols)

  # --- 3. Reformat for Galaxy-LCA Tool ---
  # The LCA tool requires specific column headers and data structure.
  # This assumes the sseqid format is: "ACCESSION###Kingdom;Phylum;...;Species"
  message("Formatting data for the LCA tool...")
  lca_formatted <- blast_out %>%
    # Split the subject ID into Accession and full Taxonomy string
    separate(sseqid, into = c("#Subject accession", "#Taxonomy"), sep = "###", remove = TRUE) %>%
    # Clean up the taxonomy string for LCA tool compatibility
    mutate(
      `#Taxonomy` = str_replace_all(`#Taxonomy`, "root_1;", ""),
      `#Taxonomy` = str_replace_all(`#Taxonomy`, ";", " / ")
    ) %>%
    # Extract the species name to create the '#Subject' column
    mutate(
      `#Subject` = str_extract(`#Taxonomy`, "(?<=/ )[^/]+$")
    ) %>%
    # Extract the taxonomy ID from the subject accession (optional, but good practice)
    mutate(
      `#Subject Taxonomy ID` = str_extract(`#Subject accession`, "^[^.]+")
    ) %>%
    # Add the database source
    mutate(
      `#Source` = db_source_name
    ) %>%
    # Rename columns to match LCA tool requirements
    rename(
      `#Query ID` = qseqid,
      `#Identity percentage` = pident,
      `#Coverage` = qcovs,
      `#evalue` = evalue,
      `#bitscore` = bitscore
    ) %>%
    # Select and reorder columns to the final required format
    select(
      `#Query ID`,
      `#Subject`,
      `#Subject accession`,
      `#Subject Taxonomy ID`,
      `#Identity percentage`,
      `#Coverage`,
      `#evalue`,
      `#bitscore`,
      `#Source`,
      `#Taxonomy`
    )

  # --- 4. Write Output and Return ---
  # Write the formatted table for the LCA tool
  message("Writing formatted LCA input file to: ", lca_input_path)
  write.table(
    lca_formatted,
    file = lca_input_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  return(lca_input_path)
}
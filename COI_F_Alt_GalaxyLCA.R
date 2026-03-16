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
  # For the new BLAST output, use the actual column names from the file
  # e.g., query_seqid, result_seqid, evalue, bitscore, length, percent_coverage, nident, percent_identity, accnumb, sciname, phylum, kingdom, class, order, family, genus, species, marker, method
  message("\nReading BLAST output from: ", blast_output_path)
  blast_out <- read.delim(blast_output_path, header = TRUE, sep = "\t", quote = "")

  # --- 3. Reformat for Galaxy-LCA Tool ---
  # The LCA tool requires specific column headers and data structure.
  # This assumes the result_seqid format is: "ACCESSION###Kingdom;Phylum;...;Species"
  message("Formatting data for the LCA tool...")
  lca_formatted <- blast_out %>%
    separate(result_seqid, into = c("#Subject accession", "#Taxonomy"), sep = "###", remove = TRUE, fill = "right") %>%
    mutate(
      `#Taxonomy` = str_replace_all(`#Taxonomy`, "root_1;", ""),
      `#Taxonomy` = str_replace_all(`#Taxonomy`, ";", " / ")
    ) %>%
    mutate(
      `#Subject` = str_extract(`#Taxonomy`, "(?<=/ )[^/]+$")
    ) %>%
    mutate(
      `#Subject Taxonomy ID` = str_extract(`#Subject accession`, "^[^.]+")
    ) %>%
    mutate(
      `#Source` = db_source_name
    ) %>%
    {
      nms <- names(.)
      nms[nms == "query_seqid"] <- "#Query ID"
      nms[nms == "percent_identity"] <- "#Identity percentage"
      nms[nms == "percent_coverage"] <- "#Coverage"
      nms[nms == "evalue"] <- "#evalue"
      nms[nms == "bitscore"] <- "#bitscore"
      names(.) <- nms
      .
    } %>%
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
    ) %>%
    mutate(
      `#Identity percentage` = as.numeric(replace_na(`#Identity percentage`, 0)),
      `#Coverage` = as.numeric(replace_na(`#Coverage`, 0)),
      `#evalue` = as.numeric(replace_na(`#evalue`, 0)),
      `#bitscore` = as.numeric(replace_na(`#bitscore`, 0))
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
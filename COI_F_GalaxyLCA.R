################################################################
# GALAXY-LCA FORMATTING - DNA metabarcoding processing
################################################################

# This function formats the BLAST output to be compatible with the
# galaxy-tool-lca python script (https://github.com/naturalis/galaxy-tool-lca).

format_blast_for_lca <- function(blast_output_path, lca_input_path, db_source_name = "MIDORI2") {

  # Load required packages
  suppressPackageStartupMessages({
    require(dplyr)
    require(tidyr)
    require(stringr)
  })

  # Define column names based on the BLAST output format
  blast_cols <- c(
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs"
  )

  # Read BLAST output
  message("Reading BLAST output from: ", blast_output_path)
  blast_out <- read.delim(blast_output_path, header = FALSE, col.names = blast_cols)

  # --- Reformat for Galaxy-LCA Tool ---
  # The LCA tool requires specific column headers and data structure.
  # We will parse the subject ID (sseqid) to extract taxonomy.
  # This assumes the sseqid format is: "ACCESSION###Kingdom;Phylum;...;Species"

  message("Formatting data for LCA tool...")
  lca_formatted <- blast_out %>%
    # Split the subject ID into Accession and full Taxonomy string
    separate(sseqid, into = c("#Subject accession", "#Taxonomy"), sep = "###", remove = TRUE) %>%
    # Clean up the taxonomy string
    mutate(
      `#Taxonomy` = str_replace_all(`#Taxonomy`, "root_1;", ""),
      `#Taxonomy` = str_replace_all(`#Taxonomy`, ";", " / ")
    ) %>%
    # Extract the species name to create the '#Subject' column
    mutate(
      `#Subject` = str_extract(`#Taxonomy`, "(?<=/ )[^/]+$")
    ) %>%
    # Extract the taxonomy ID from the subject accession
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
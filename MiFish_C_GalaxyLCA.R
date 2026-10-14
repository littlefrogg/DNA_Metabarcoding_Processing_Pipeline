################################################################
# (MiFish-C) PREPARE FOR GALAXY-LCA - DNA metabarcoding processing
################################################################
# Purpose: Format BLAST output to be compatible with the galaxy-tool-lca
#          python script for Lowest Common Ancestor (LCA) taxonomic assignment.
#
# This version is tolerant of two BLAST output styles:
# - headered TSV (column names present, e.g., produced by ALT BLAST)
# - headerless TSV using the COI/standard outfmt (reads into defined columns)
# It also expands multi-hit `result_seqid` cells separated by '-' into separate rows.
################################################################

format_blast_for_lca <- function(blast_output_path, lca_input_path, db_source_name) {

  suppressPackageStartupMessages({
    require(dplyr)
    require(tidyr)
    require(stringr)
  })

  if (!file.exists(blast_output_path) || file.info(blast_output_path)$size == 0) {
    stop("BLAST output file not found or is empty at: ", blast_output_path)
  }

  message("\nReading BLAST output from: ", blast_output_path)
  # Peek at first line to detect header presence
  first_line <- readLines(blast_output_path, n = 1)
  has_header <- grepl("qseqid|query_seqid|query_seqid", first_line, ignore.case = TRUE)

  if (has_header) {
    blast_out <- read.delim(blast_output_path, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  } else {
    blast_cols <- c(
      "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
      "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs"
    )
    blast_out <- read.delim(blast_output_path, header = FALSE, sep = "\t", quote = "", col.names = blast_cols, stringsAsFactors = FALSE)
  }

  # Normalize names
  if ("sseqid" %in% names(blast_out)) names(blast_out)[names(blast_out) == "sseqid"] <- "result_seqid"
  if ("qseqid" %in% names(blast_out)) names(blast_out)[names(blast_out) == "qseqid"] <- "query_seqid"
  if ("pident" %in% names(blast_out)) names(blast_out)[names(blast_out) == "pident"] <- "percent_identity"
  if ("qcovs" %in% names(blast_out)) names(blast_out)[names(blast_out) == "qcovs"] <- "percent_coverage"

  # Expand multi-hit result_seqid cells separated by '-'
  if ("result_seqid" %in% names(blast_out)) {
    blast_out <- blast_out %>%
      mutate(result_seqid = as.character(result_seqid)) %>%
      tidyr::separate_rows(result_seqid, sep = "-")
  }

  message("Formatting data for the LCA tool...")
  lca_formatted <- blast_out %>%
    tidyr::separate(result_seqid, into = c("#Subject accession", "#Taxonomy"), sep = "###", remove = TRUE, fill = "right") %>%
    mutate(
      `#Taxonomy` = str_replace_all(`#Taxonomy`, "root_1;", ""),
      `#Taxonomy` = str_replace_all(`#Taxonomy`, ";", " / ")
    ) %>%
    mutate(
      `#Subject` = str_extract(`#Taxonomy`, "(?<=/ )[^/]+$"),
      `#Subject Taxonomy ID` = str_extract(`#Subject accession`, "^[^.]+"),
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

  message("Writing formatted LCA input file to: ", lca_input_path)
  write.table(lca_formatted, file = lca_input_path, sep = "\t", row.names = FALSE, quote = FALSE)
  return(lca_input_path)
}

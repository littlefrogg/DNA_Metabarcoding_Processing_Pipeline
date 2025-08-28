#Step D script for eDNA sequence analysis pipeline #
# Paige Smallman, 2025
#Creating phylosec object

create_phyloseq <- function(
    otu_table_rds,
    tax_table_rds,
    metadata_path,
    output_dir,
    project_id,
    sample_id_col = 2,
    fix_mismatches = FALSE) {
  
# Load required packages
  suppressPackageStartupMessages({
    require(phyloseq)
    require(here)
    require(tidyverse)
  })
  
# INPUT VALIDATION

# Check input files exist
inputs <- c(otu_table_rds, tax_table_rds, metadata_path)
  if(any(!file.exists(inputs))) {
    missing <- inputs[!file.exists(inputs)]
    stop("Missing input files:\n", paste(missing, collapse = "\n"))
  }
  
# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
# DATA LOADING
  
message("\nLoading input data...")
  
# Load OTU table
otu_table <- readRDS(otu_table_rds)
message("OTU table dimensions: ", paste(dim(otu_table), collapse = " x "))
  
# Load taxonomy table
tax_table <- readRDS(tax_table_rds)
message("Taxonomy table dimensions: ", paste(dim(tax_table), collapse = " x "))
  
# Validate table compatibility
  if(nrow(otu_table) != nrow(tax_table)) {
    stop("OTU/Taxonomy table mismatch: ",
         nrow(otu_table), " OTUs vs ",
         nrow(tax_table), " taxa")
  }
  
# Load metadata
metadata <- read.csv(metadata_path, header = TRUE, sep = ",", row.names = sample_id_col)
message("Metadata dimensions: ", paste(dim(metadata), collapse = " x "))
  
# SAMPLE VALIDATION
  
# Get sample IDs
otu_samples <- rownames(otu_table)
meta_samples <- rownames(metadata)
  
  mismatches <- list(
    missing_in_meta = setdiff(otu_samples, meta_samples),
    missing_in_otu = setdiff(meta_samples, otu_samples)
  )
  
if(length(unlist(mismatches)) > 0) {
    message("\nSample mismatches detected:")
    message("- Samples in OTU table missing from metadata: ", 
            length(mismatches$missing_in_meta))
    message("- Samples in metadata missing from OTU table: ", 
            length(mismatches$missing_in_otu))
    
    if(fix_mismatches) {
      message("\nResolving mismatches...")
      # Subset to intersecting samples
      common_samples <- intersect(otu_samples, meta_samples)
      
      otu_table <- otu_table[ ,common_samples]
      metadata <- metadata[common_samples, ]
      
      message("Retained ", length(common_samples), " matching samples")
    } else {
      warning("Sample mismatches present. Use fix_mismatches=TRUE to automatically subset")
    }
  }
  
# PHYLOSEQ CONSTRUCTION
  
message("\nConstructing phyloseq object...")
  
  ps <- phyloseq(
    otu_table(otu_table, taxa_are_rows = TRUE),
    sample_data(metadata),
    tax_table(tax_table)
  )
  
message("\nPhyloseq object created:")
print(ps)
  
# OUTPUT GENERATION
  
output_path <- file.path(output_dir, paste0(project_id, "_phyloseq.rds"))
saveRDS(ps, output_path)
message("\nSaved phyloseq object to: ", output_path)
  
# DIAGNOSTIC OUTPUTS
  
# Save sample data for verification
  write.csv(
    as(sample_data(ps), "data.frame"),
    file = file.path(output_dir, paste0(project_id, "_phy_samples.csv")),
    row.names = TRUE
  )
  
  return(list(
    phyloseq = ps,
    output_path = output_path,
    mismatches = mismatches
  ))
}

# Example usage:
# phy_results <- create_phyloseq(
#   otu_table_rds = here("outputs/ROHR05/ROHR05.nochim.rds"),
#   tax_table_rds = here("outputs/ROHR05/ROHR05_taxtable.rds"),
#   metadata_path = here("inputs/metadata/metadata.csv"),
#   output_dir = here("outputs/ROHR05"),
#   project_id = "ROHR05",
#   sample_id_col = 2,
#   fix_mismatches = TRUE
# )
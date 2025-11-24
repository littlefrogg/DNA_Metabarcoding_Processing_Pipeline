################################################################
# CLUSTER ASVs COI - DNA metabarcoding processing
################################################################

# Cluster COI ASVs using DECIPHER
cluster_coi_asvs <- function(ps_object, 
                             output_dir = output_root,
                             project_id = project_id,
                             cutoff = 0.03,         # 97% similarity
                             min_coverage = 0.8,
                             processors = NULL) {
  
# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
# Extract sequences from phyloseq
dna <- refseq(ps_object)
  if(is.null(dna)) stop("Phyloseq object missing reference sequences")
  
# Cluster ASVs
message("Found ", length(dna), " ASVs.")
message("\n[", Sys.time(), "] Starting ASV clustering...")
  clusters <- Clusterize(
    dna,
    cutoff = cutoff,
    minCoverage = min_coverage,
    processors = processors,
    verbose = TRUE
  )
message("Clustering complete. Resulted in ", nrow(clusters), " OTUs.")

# Process cluster assignments
  cluster_df <- tibble(
    ASV = names(dna),
    OTU = paste0("OTU_", clusters$cluster)
  )
  
# Save cluster assignments
cluster_file <- file.path(output_dir, paste0(project_id, "_clusters.csv"))
write_csv(cluster_df, cluster_file)
message("Saved cluster assignments to: ", cluster_file)
  
# Create OTU table
message("\nCreating OTU table...")
otu <- as(otu_table(ps_object), "matrix")
if (!taxa_are_rows(ps_object)) otu <- t(otu)
cluster_df <- cluster_df[match(rownames(otu), cluster_df$ASV), ]
otu_mat <- rowsum(otu, group = cluster_df$OTU)
otu_mat <- t(otu_mat) # if you want OTUs as columns
stopifnot(all(rownames(otu) == cluster_df$ASV))

# Create consensus taxonomy
message("Generating consensus taxonomy...")
tax_df <- tax_table(ps_object) %>%
  as.data.frame() %>%
  rownames_to_column("ASV") %>%
  left_join(cluster_df, by = "ASV") %>%
  group_by(OTU) %>%
  summarize(across(where(is.character), ~{
    vals <- na.omit(.)
    if (length(vals) == 0) NA_character_ else names(which.max(table(vals)))
  })) %>%
  column_to_rownames("OTU") %>%
  as.matrix()
  
# Build new phyloseq object
message("Constructing phyloseq object...")
  ps_otu <- phyloseq(
    otu_table(otu_mat, taxa_are_rows = FALSE),
    tax_table(tax_df),
    sample_data(ps_object)
  )
  
# Save results
output_file <- file.path(output_dir, paste0(project_id, "_clustered_phyloseq.rds"))
saveRDS(ps_otu, output_file)
message("\nSuccessfully saved clustered phyloseq to: ", output_file)
  
return(ps_otu)
}


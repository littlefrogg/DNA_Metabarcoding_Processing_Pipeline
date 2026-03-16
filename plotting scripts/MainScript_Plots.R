################################################################
# PLOTS MAIN SCRIPT - DNA metabarcoding visualization
################################################################

# COI Plots
# Jan 30, 2025
# R version 4.5

################################################################
# (1) PACKAGE MANAGEMENT
################################################################

# --- Set Working Directory ---
#setwd("/home/paiges/scratch/eDNA_COI/eDNA")

# Load necessary libraries
cran_packages <- c(
  "here", "ggplot2", "dplyr", "remotes", "tidyverse"
)
bioc_packages <- c(
  "phyloseq"
)
github_packages <- list("pairwiseAdonis")

# Install CRAN packages
to_install_cran <- setdiff(cran_packages, rownames(installed.packages()))
if(length(to_install_cran)) install.packages(to_install_cran)

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- setdiff(bioc_packages, rownames(installed.packages()))
if(length(to_install_bioc)) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

# Install GitHub packages if not present
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
for(pkg in names(github_packages)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    remotes::install_github(github_packages[[pkg]])
  }
}

# Load all packages
all_packages <- unique(c(cran_packages, bioc_packages, names(github_packages)))
invisible(lapply(all_packages, function(pkg) library(pkg, character.only = TRUE)))

################################################################
# (2) USER CONFIGURATION - SET ALL PARAMETERS HERE
################################################################

# --- Project & File Settings ---
#project_ids   <- c("COI_plate4", "COI_plate3", "COI_plate2", "COI_plate1", "COI_2020_21") # Vector of project IDs to combine (e.g. c("Plate1", "Plate2"))
#run_name      <- "ALL_COI"    # Name for this analysis run (folder name for figures)
#project_ids   <- c("COI_plate4", "COI_plate3", "COI_plate2", "COI_plate1") # Vector of project IDs to combine (e.g. c("Plate1", "Plate2"))
#run_name      <- "CombinedCOI"    # Name for this analysis run (folder name for figures)
project_ids   <- c("COI_plate4_alt", "COI_plate3_alt", "COI_plate2_alt", "COI_plate1_alt", "COI_2020_21_alt") # Vector of project IDs to combine (e.g. c("Plate1", "Plate2"))
run_name      <- "ALL_COI_alt"    # Name for this analysis run (folder name for figures)
primer        <- "COI"
metadata_file <- "metadata_unique.csv"
tax_db        <- "MIDORI2_UNIQ_NUC_GB267_CO1_DADA2.fasta"

################################################################
# (3) AUTOMATED SETUP - DO NOT EDIT BELOW THIS LINE
################################################################

# --- Define Directory and File Paths ---
# Core directories
input_root   <- here("inputs")
scripts_root <- here("scripts")
output_root  <- here("outputs", project_ids)
figures_root <- here("figures", run_name)
tax_root     <- here("tax", tax_db)

# Input paths
metadata_path <- here(input_root, "metadata", metadata_file)

# Merging Phyloseq objects
if (!exists("ps_final")) {
  message("Loading final phyloseq objects for projects: ", paste(project_ids, collapse = ", "))
  ps_list <- list()
  for (pid in project_ids) {
    # Construct path for this project
    fpath <- here("outputs", pid, paste0(pid, "_final_phyloseq.rds"))
    if (file.exists(fpath)) {
      message("  Loading: ", fpath)
      ps_list[[pid]] <- readRDS(fpath)
    } else {
      warning("  File not found: ", fpath)
    }
  }
  if (length(ps_list) == 0) {
    stop("No phyloseq objects loaded. Check project IDs and paths.")
  } else if (length(ps_list) == 1) {
    ps_final <- ps_list[[1]]
  } else {
    message("Merging ", length(ps_list), " phyloseq objects...")
    # Iteratively merge
    ps_final <- ps_list[[1]]
    for (i in 2:length(ps_list)) {
      ps_final <- merge_phyloseq(ps_final, ps_list[[i]])
    }
  }
}

# Remove empty samples (all-zero counts) to create ps_final_no_empty
ps_final_no_empty <- prune_samples(sample_sums(ps_final) > 0, ps_final)
if (any(sample_sums(ps_final) == 0)) {
  message("Removed ", sum(sample_sums(ps_final) == 0), " empty samples from ps_final. Result stored in ps_final_no_empty.")
} else {
  message("No empty samples detected in ps_final.")
}

# Clean up tax table for ps_final_no_empty (remove trailing numbers)
if ("tax_table" %in% slotNames(ps_final_no_empty)) {
  taxmat_no_empty <- as.matrix(tax_table(ps_final_no_empty))
  taxmat_no_empty_clean <- apply(taxmat_no_empty, 2, function(x) gsub("_[0-9]+$", "", x))
  tax_table(ps_final_no_empty) <- taxmat_no_empty_clean
  message("Cleaned trailing numbers from tax_table in ps_final_no_empty.")
}

# Clean up tax table

# Remove trailing numbers (e.g., _6340) from all taxonomic ranks in tax_table
if ("tax_table" %in% slotNames(ps_final)) {
  taxmat <- as.matrix(tax_table(ps_final))
  taxmat_clean <- apply(taxmat, 2, function(x) gsub("_[0-9]+$", "", x))
  tax_table(ps_final) <- taxmat_clean
  message("Cleaned trailing numbers from tax_table.")
}

################################################################
# (3b) DATA INTEGRITY CHECKS BEFORE PLOTTING
################################################################

# Check for samples with NA in key metadata columns
meta <- as.data.frame(sample_data(ps_final))
key_cols <- c("Year", "season") # Adjust as needed for your grouping variables
for (col in key_cols) {
  if (col %in% colnames(meta)) {
    n_na <- sum(is.na(meta[[col]]))
    if (n_na > 0) {
      warning(sprintf("%d samples have NA in '%s'. These may appear as NA categories in plots.", n_na, col))
    }
  } else {
    warning(sprintf("Column '%s' not found in sample_data. Check your metadata.", col))
  }
}

# Check for samples in phyloseq not in metadata
meta_samples <- rownames(meta)
phylo_samples <- sample_names(ps_final)
missing_in_meta <- setdiff(phylo_samples, meta_samples)
if (length(missing_in_meta) > 0) {
  warning(sprintf("%d samples in phyloseq object not found in metadata: %s", length(missing_in_meta), paste(missing_in_meta, collapse=", ")))
}

# Check for all-zero samples
zero_samples <- sample_names(ps_final)[sample_sums(ps_final) == 0]
if (length(zero_samples) > 0) {
  warning(sprintf("%d samples have all zero counts and may cause NA categories in plots: %s", length(zero_samples), paste(zero_samples, collapse=", ")))
}

# Drop unused factor levels in metadata
for (col in key_cols) {
  if (col %in% colnames(meta) && is.factor(meta[[col]])) {
    meta[[col]] <- droplevels(meta[[col]])
  }
}
sample_data(ps_final) <- meta

################################################################
# 4 DIVERSITY SUMMARY PLOTS
################################################################

# --- Define Palette ---
Phylum_palette <- c(
  "Porifera" = "#EEC644", "Cnidaria" = "#FF7F50", 
  "Mollusca" = "#59d99dff", "Annelida" = "#54D164", "Platyhelminthes" = "#7CD178", 
  "Arthropoda" = "#4596D9", "Gastrotricha" = "#7764D6", "Nemertea" = "#DC6160", 
  "Echinodermata" = "#B8426F",
  "Bryozoa" = "#A9DC6A",
  "Chaetognatha" = "#5E7DD8",
  "Ctenophora" = "#E7635E",
  "Hemichordata" = "#7F6B9E",
  "Nematoda" = "#A75B89",
  "Xenacoelomorpha" = "#8C5A8F",
  # Invertebrate chordates (tunicates: Ascidiacea, Thaliacea)
  "Chordata" = "#9B59B6"
)

# Centralized OTU color styling (edit here)
otu_color_config <- list(
  use_base_for_phylum = TRUE,      # keep exact base phylum colors
  lighten_non_phylum = 0.18,       # lighten Class/Order/Family/Genus/Species bars
  bar_alpha_phylum = 1.00,
  bar_alpha_other = 0.85,
  bar_alpha_highest_rank = 0.88,
  highest_rank_light_low = 0.45,
  highest_rank_light_high = 0.15,
  fallback_color = "#BDBDBD"
)

# --- Debugging: Check Phylum levels and palette matching ---
taxmat_plot <- as.data.frame(tax_table(ps_final_no_empty))
if ("Phylum" %in% colnames(taxmat_plot)) {
  phyla_present <- sort(unique(as.character(taxmat_plot$Phylum)))
  message("Phyla present in ps_final_no_empty: ", paste(phyla_present, collapse=", "))
  palette_names <- names(Phylum_palette)
  missing_in_palette <- setdiff(phyla_present, palette_names)
  missing_in_data <- setdiff(palette_names, phyla_present)
  if (length(missing_in_palette) > 0) {
    warning(sprintf("Phyla in data not in palette: %s", paste(missing_in_palette, collapse=", ")))
  }
  if (length(missing_in_data) > 0) {
    message(sprintf("Palette contains colors for phyla not present in data: %s", paste(missing_in_data, collapse=", ")))
  }
  # Optionally, restrict palette to only present phyla for plotting
  Phylum_palette_plot <- Phylum_palette[names(Phylum_palette) %in% phyla_present]
} else {
  warning("Phylum column not found in tax_table. Barplot colors may not work.")
  Phylum_palette_plot <- Phylum_palette
}

# --- Load Helper Functions ---
source(file.path(scripts_root, "taxa_color_utils.R"))

# --- Stats ---

# Permanova
source(file.path(scripts_root, "COI_H_PERMANOVA.R"))

# --- Subset to Invertebrates for PERMANOVA ---
message("Subsetting for Invertebrate diversity for PERMANOVA...")

# Subset to invertebrate phyla for plotting (only those in palette, after cleaning)
# Include tunicates (Ascidiacea, Thaliacea) which are invertebrate chordates
ps_inv_no_empty <- subset_taxa(
  ps_final_no_empty,
  !is.na(Phylum) &
  (
    !(Phylum %in% c("Chordata", "Vertebrata")) |  # Non-chordates
    (Phylum == "Chordata" & Class %in% c("Ascidiacea", "Thaliacea"))  # Invertebrate chordates (tunicates)
  )
)

# Clean up tax table for ps_inv_no_empty (remove trailing numbers)
if ("tax_table" %in% slotNames(ps_inv_no_empty)) {
  taxmat_inv_no_empty <- as.matrix(tax_table(ps_inv_no_empty))
  taxmat_inv_no_empty_clean <- apply(taxmat_inv_no_empty, 2, function(x) gsub("_[0-9]+$", "", x))
  tax_table(ps_inv_no_empty) <- taxmat_inv_no_empty_clean
  message("Cleaned trailing numbers from tax_table in ps_inv_no_empty.")
  # Debug: Check Phylum column after cleaning
  cat("[DEBUG] Phylum column after cleaning trailing numbers (head):\n")
  print(head(as.data.frame(tax_table(ps_inv_no_empty))$Phylum))
}

# Filter to only taxa with Phylum in palette
taxmat_inv <- as.data.frame(tax_table(ps_inv_no_empty))
if ("Phylum" %in% colnames(taxmat_inv)) {
  in_palette <- taxmat_inv$Phylum %in% names(Phylum_palette)
  keep_taxa <- rownames(taxmat_inv)[in_palette]
  ps_inv_no_empty <- prune_taxa(keep_taxa, ps_inv_no_empty)
  taxmat_inv <- as.data.frame(tax_table(ps_inv_no_empty))
  # Assign factor levels for Phylum (for plotting)
  invertebrate_phyla <- names(Phylum_palette)
  taxmat_inv$Phylum <- factor(as.character(taxmat_inv$Phylum), levels = invertebrate_phyla)
  # IMPORTANT: Convert factor back to character before assigning to tax_table (matrix)
  tax_table(ps_inv_no_empty)[, "Phylum"] <- as.character(taxmat_inv$Phylum)
  Phylum_palette_inv <- Phylum_palette
  message("Invertebrate phyla (palette): ", paste(invertebrate_phyla, collapse=", "))
  message("Invertebrate phyla present: ", paste(sort(unique(as.character(taxmat_inv$Phylum))), collapse=", "))
  # Debug: Check Phylum column after factor assignment
  cat("[DEBUG] Phylum column after factor assignment (head):\n")
  print(head(as.data.frame(tax_table(ps_inv_no_empty))$Phylum))
  cat("[DEBUG] Unique phyla after filtering and factor assignment:\n")
  print(sort(unique(as.character(taxmat_inv$Phylum))))
  cat("[DEBUG] NAs in Phylum column after filtering:", sum(is.na(taxmat_inv$Phylum)), "\n")
  cat("[DEBUG] Empty strings in Phylum column after filtering:", sum(taxmat_inv$Phylum == "", na.rm=TRUE), "\n")
} else {
  warning("Phylum column not found in tax_table for invertebrates.")
  Phylum_palette_inv <- Phylum_palette
}

permanova_res <- run_permanova_richness(
  ps_inv_no_empty,
  nperm = 2999,
  seed = 123,
  results_file = file.path(figures_root, "permanova_results.txt")
)

# Extract p-value and R² from the PERMANOVA table
permanova_p <- NA
permanova_r2 <- NA
if (inherits(permanova_res$permanova, "anova")) {
  # Try to get values from the table (handles both with and without "Model" row)
  if (nrow(permanova_res$permanova) > 0) {
    # If there's a "Model" row, use it; otherwise use the first term row
    if ("Model" %in% rownames(permanova_res$permanova)) {
      permanova_p <- permanova_res$permanova["Model", "Pr(>F)"]
      permanova_r2 <- permanova_res$permanova["Model", "R2"]
    } else if (!is.na(permanova_res$permanova[1, "Pr(>F)"]) && !is.na(permanova_res$permanova[1, "R2"])) {
      # Use first row if available
      permanova_p <- permanova_res$permanova[1, "Pr(>F)"]
      permanova_r2 <- permanova_res$permanova[1, "R2"]
    }
  }
}

# --- Plots ---

# Summary Plots (Phylum Richness, Alpha Diversity)
source(file.path(scripts_root, "COI_H_SummaryPlots.R"))
generate_summary_plots(
  ps_object = ps_final_no_empty,
  figures_dir = figures_root,
  project_id = run_name,
  tax_palette = Phylum_palette
)

# Taxa Plots
source(file.path(scripts_root, "COI_H_TaxaPlots.R"))
generate_taxa_plots(
  ps_object   = ps_final_no_empty,
  figures_dir = figures_root,
  project_id  = run_name,
  tax_palette = Phylum_palette
)

# PCoA Plots
source(file.path(scripts_root, "COI_H_PCOA_plots.R"))
# Eukaryote PCoA
plot_pcoa_jaccard(
  ps_object = ps_final_no_empty,
  figures_dir = figures_root,
  project_id = run_name
)
# Invertebrate PCoA
plot_pcoa_jaccard(
  ps_object = ps_inv_no_empty,
  figures_dir = file.path(figures_root),
  project_id = paste0(run_name, "_invertebrates")
)

# Boxplots for alpha diversity
source(file.path(scripts_root, "COI_H_Boxplots.R"))

# --- Subset to Invertebrates for Boxplots ---
message("Subsetting for Invertebrate diversity for boxplots...")
# Include tunicates (Ascidiacea, Thaliacea) which are invertebrate chordates
ps_inv_no_empty_box <- subset_taxa(
  ps_final_no_empty,
  !is.na(Phylum) &
  (
    !(Phylum %in% c("Chordata", "Vertebrata")) |  # Non-chordates
    (Phylum == "Chordata" & Class %in% c("Ascidiacea", "Thaliacea"))  # Invertebrate chordates (tunicates)
  )
)

# Clean up tax table for ps_inv_no_empty_box (remove trailing numbers)
if ("tax_table" %in% slotNames(ps_inv_no_empty_box)) {
  taxmat_inv_no_empty_box <- as.matrix(tax_table(ps_inv_no_empty_box))
  taxmat_inv_no_empty_box_clean <- apply(taxmat_inv_no_empty_box, 2, function(x) gsub("_[0-9]+$", "", x))
  tax_table(ps_inv_no_empty_box) <- taxmat_inv_no_empty_box_clean
  message("Cleaned trailing numbers from tax_table in ps_inv_no_empty_box.")
}

# Ensure Phylum is a factor with only present levels for invertebrates (boxplots)
taxmat_inv_box <- as.data.frame(tax_table(ps_inv_no_empty_box))
if ("Phylum" %in% colnames(taxmat_inv_box)) {
  invertebrate_phyla <- names(Phylum_palette)
  taxmat_inv_box$Phylum <- factor(as.character(taxmat_inv_box$Phylum), levels = invertebrate_phyla)
  # IMPORTANT: Convert factor back to character before assigning to tax_table (matrix)
  tax_table(ps_inv_no_empty_box)[, "Phylum"] <- as.character(taxmat_inv_box$Phylum)
  Phylum_palette_inv_box <- Phylum_palette
  message("Invertebrate phyla (palette, boxplots): ", paste(invertebrate_phyla, collapse=", "))
  message("Invertebrate phyla present (boxplots): ", paste(sort(unique(as.character(taxmat_inv_box$Phylum))), collapse=", "))
} else {
  warning("Phylum column not found in tax_table for invertebrates (boxplots).")
  Phylum_palette_inv_box <- Phylum_palette
}

generate_alpha_boxplots(
  ps_object = ps_inv_no_empty_box,
  figures_dir = figures_root,
  project_id = run_name,
  permanova_p = permanova_p,
  permanova_r2 = permanova_r2,
  pairwise_results = permanova_res$pairwise
)

# OTU Barplots (labeled with highest taxonomy)
source(file.path(scripts_root, "COI_H_OTU_Barplots.R"))
taxmat_otu <- as.data.frame(tax_table(ps_inv_no_empty))

# OTU barplots for invertebrates only, using invertebrate palette
generate_otu_barplots(
  ps_object   = ps_inv_no_empty,
  figures_dir = figures_root,
  project_id  = paste0(run_name, "_invertebrates"),
  tax_palette = Phylum_palette_inv,
  color_config = otu_color_config
)

# Rarefaction Curves
source(file.path(scripts_root, "COI_H_Rarefaction_plots.R"))
plot_rarefaction_curves(
  ps_object = ps_final_no_empty,
  figures_dir = figures_root,
  project_id = run_name
)

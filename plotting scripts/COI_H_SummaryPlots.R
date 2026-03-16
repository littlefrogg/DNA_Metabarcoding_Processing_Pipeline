################################################################
# (COI-H) SUMMARY PLOTS - DNA metabarcoding processing
################################################################
# Purpose: Generate quick summary plots for the final phyloseq object.
#
# Input:
#   - ps_object:    The final phyloseq object.
#   - figures_dir:  Directory to save the plots.
#   - project_id:   Project identifier for file naming.
#
# Output:
#   - PNG files of taxonomic composition, alpha diversity, and beta diversity.
#
# Author: Paige Smallman, 2025
################################################################
generate_summary_plots <- function(ps_object, figures_dir, project_id, tax_palette) {
  # --- 1. Load Packages ---
  suppressPackageStartupMessages({
    require(phyloseq)
    require(ggplot2)
    require(dplyr)
    require(tibble)
    require(colorspace) # Required for lightening colors
  })
  
  clean_grouping <- function(x) {
    x <- trimws(as.character(x))
    x <- tolower(x)
    x <- gsub("_", " ", x)
    x <- tools::toTitleCase(x)
    x
  }
  
  # --- 2. Data Prep ---
  # Clean grouping columns in sample_data
  for (col in c("season", "local", "site", "sample_type")) {
    if (col %in% colnames(sample_data(ps_object))) {
      sample_data(ps_object)[[col]] <- clean_grouping(sample_data(ps_object)[[col]])
    }
  }
  
  if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
  message("\nGenerating summary plots in: ", figures_dir)
  
  # --- 2a. Taxonomic Composition ---
  # Generic function to plot richness by any taxonomic level and grouping variable
  plot_taxa_richness <- function(ps, group_var, tax_level = "Phylum", filename, title_ext) {
    message("- Generating ", tax_level, " richness barplot by ", group_var, "...")
    ps_merged <- merge_samples(ps, group_var)
    sample_data(ps_merged)[[group_var]] <- rownames(sample_data(ps_merged))
    otu_table(ps_merged) <- ifelse(otu_table(ps_merged) > 0, 1, 0)
    
    # Get consistent colors
    taxa_cols <- get_taxa_colors(ps_merged, level = tax_level, base_palette = tax_palette)
    
    p <- plot_bar(ps_merged, x = group_var, fill = tax_level) +
      geom_bar(aes_string(color = tax_level, fill = tax_level), stat = "identity", position = "stack") +
      scale_fill_manual(values = taxa_cols) +
      scale_color_manual(values = taxa_cols) +
      theme_minimal() +
      labs(title = paste(project_id, "-", title_ext), y = "Richness (Number of OTUs)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(figures_dir, filename), p, width = 10, height = 7)
  }
  
  # Phylum-level richness plots
  plot_taxa_richness(ps_object, "site", "Phylum", paste0(project_id, "_taxa_by_site_richness.png"), "Phylum Richness by Site")
  plot_taxa_richness(ps_object, "season", "Phylum", paste0(project_id, "_taxa_by_season_richness.png"), "Phylum Richness by Season")
  plot_taxa_richness(ps_object, "local", "Phylum", paste0(project_id, "_taxa_by_location_richness.png"), "Phylum Richness by Location")
  
  # Class-level richness plots
  plot_taxa_richness(ps_object, "site", "Class", paste0(project_id, "_class_by_site_richness.png"), "Class Richness by Site")
  plot_taxa_richness(ps_object, "season", "Class", paste0(project_id, "_class_by_season_richness.png"), "Class Richness by Season")
  plot_taxa_richness(ps_object, "local", "Class", paste0(project_id, "_class_by_location_richness.png"), "Class Richness by Location")
  
  # --- Alpha Diversity ---
  # Note: Alpha diversity uses sample metadata colors, not taxonomic colors
  message("- Generating alpha diversity plots...")
  p_alpha_local <- plot_richness(ps_object, x = "local", measures = "Observed") +
    geom_boxplot(aes(fill = local), alpha = 0.7) +
    theme_bw() + 
    labs(title = "Alpha Diversity",  y = "Richness (Number of OTUs)")
  ggsave(file.path(figures_dir, paste0(project_id, "_alpha_location.png")), p_alpha_local, width = 8, height = 6)

  # Note: Alpha diversity uses sample metadata colors, not taxonomic colors
  p_alpha_season <- plot_richness(ps_object, x = "season", measures = "Observed") +
    geom_boxplot(aes(fill = season), alpha = 0.7) +
    theme_bw() + 
    labs(title = "Alpha Diversity",  y = "Richness (Number of OTUs)")
  ggsave(file.path(figures_dir, paste0(project_id, "_alpha_season.png")), p_alpha_season, width = 8, height = 6)
  
  message("Summary plots generated successfully.")
}
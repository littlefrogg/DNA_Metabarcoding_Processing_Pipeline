################################################################
# (COI-H) PCOA PLOTS - DNA metabarcoding detailed visualization
################################################################
# Purpose: Generate PCoA plots for the final phyloseq object.
#
# Author: Paige Smallman, 2026
################################################################

plot_pcoa_jaccard <- function(ps_object, figures_dir, project_id) {
  suppressPackageStartupMessages({
    require(phyloseq)
    require(ggplot2)
  })

  clean_grouping <- function(x) {
    x <- trimws(as.character(x))
    x <- tolower(x)
    x <- gsub("_", " ", x)
    x <- tools::toTitleCase(x)
    x
  }

  # Remove samples with zero taxa (empty rows cause NA distances)
  sample_sums_vec <- sample_sums(ps_object)
  empty_samples <- names(sample_sums_vec[sample_sums_vec == 0])
  
  if (length(empty_samples) > 0) {
    message("Removing ", length(empty_samples), " empty samples with zero reads: ", 
            paste(head(empty_samples, 5), collapse = ", "),
            ifelse(length(empty_samples) > 5, "...", ""))
    ps_object <- prune_samples(sample_sums(ps_object) > 0, ps_object)
  }
  
  # Also remove taxa with zero counts across remaining samples
  ps_object <- prune_taxa(taxa_sums(ps_object) > 0, ps_object)
  
  message("After filtering: ", nsamples(ps_object), " samples, ", ntaxa(ps_object), " taxa")
  
  # Clean grouping columns in sample_data
  for (col in c("season", "local", "site", "sample_type")) {
    if (col %in% colnames(sample_data(ps_object))) {
      sample_data(ps_object)[[col]] <- clean_grouping(sample_data(ps_object)[[col]])
    }
  }
  
  # Helper function to add ellipses only for groups with enough samples
  add_ellipse_if_valid <- function(plot_obj, grouping_var, ps_obj, min_samples = 4) {
    group_counts <- table(sample_data(ps_obj)[[grouping_var]])
    valid_groups <- names(group_counts[group_counts >= min_samples])
    
    if (length(valid_groups) > 0) {
      # Filter data to only include valid groups for ellipse
      plot_obj <- plot_obj + 
        stat_ellipse(data = ~ subset(.x, .x[[grouping_var]] %in% valid_groups),
                     na.rm = TRUE)
    }
    return(plot_obj)
  }

  message("Generating sample type PCoA plot (Jaccard) in: ", figures_dir)
  tryCatch({
    ord_pcoa <- ordinate(ps_object, method = "PCoA", distance = "jaccard", binary = TRUE)
    p_beta <- plot_ordination(ps_object, ord_pcoa, color = "sample_type") +
      geom_point(size = 4) +
      theme_minimal() +
      labs(title = paste(project_id, "- PCoA (Jaccard)"))
    
    # Only add ellipses for groups with >= 4 samples
    group_counts <- table(sample_data(ps_object)$sample_type)
    if (any(group_counts >= 4)) {
      p_beta <- p_beta + stat_ellipse(na.rm = TRUE)
    }
    
    ggsave(file.path(figures_dir, paste0(project_id, "_sampletype_pcoa_jaccard.png")), p_beta, width = 8, height = 6)
  }, error = function(e) {
    message("Warning: Could not generate Beta Diversity plot: ", e$message)
  })

message("Generating season PCoA plot (Jaccard) in: ", figures_dir)
  tryCatch({
    ord_pcoa <- ordinate(ps_object, method = "PCoA", distance = "jaccard", binary = TRUE)
    p_beta <- plot_ordination(ps_object, ord_pcoa, color = "season") +
      geom_point(size = 4) +
      theme_minimal() +
      labs(title = paste(project_id, "- PCoA (Jaccard)"))
    
    group_counts <- table(sample_data(ps_object)$season)
    if (any(group_counts >= 4)) {
      p_beta <- p_beta + stat_ellipse(na.rm = TRUE)
    }
    
    ggsave(file.path(figures_dir, paste0(project_id, "_season_pcoa_jaccard.png")), p_beta, width = 8, height = 6)
  }, error = function(e) {
    message("Warning: Could not generate Beta Diversity plot: ", e$message)
  })

  message("Generating Location PCoA plot (Jaccard) in: ", figures_dir)
  tryCatch({
    ord_pcoa <- ordinate(ps_object, method = "PCoA", distance = "jaccard", binary = TRUE)
    p_beta <- plot_ordination(ps_object, ord_pcoa, color = "local") +
      geom_point(size = 4) +
      theme_minimal() +
      labs(title = paste(project_id, "- PCoA (Jaccard)"))
    
    group_counts <- table(sample_data(ps_object)$local)
    if (any(group_counts >= 4)) {
      p_beta <- p_beta + stat_ellipse(na.rm = TRUE)
    }
    
    ggsave(file.path(figures_dir, paste0(project_id, "_location_pcoa_jaccard.png")), p_beta, width = 8, height = 6)
  }, error = function(e) {
    message("Warning: Could not generate Beta Diversity plot: ", e$message)
  })

  message("Generating PCoA (Jaccard): Shapes = Location, Colors = Year")
  tryCatch({
    # 1. Ensure Year and local are factors for discrete plotting
    # We use the clean_grouping helper for the location labels
    sample_data(ps_object)$local <- factor(clean_grouping(sample_data(ps_object)$local))
    sample_data(ps_object)$Year  <- factor(sample_data(ps_object)$Year)
    
    # 2. Ordination
    ord_pcoa <- ordinate(ps_object, method = "PCoA", distance = "jaccard", binary = TRUE)
    
    # 3. Plot
    p_beta_styled <- plot_ordination(ps_object, ord_pcoa, 
                                     color = "Year", 
                                     shape = "local") +
      geom_point(size = 5, alpha = 0.8) +
      theme_minimal() +
      # stat_ellipse(aes(group = local), linetype = 2, color = "gray60") + # Optional: Group by site
      # stat_ellipse(aes(group = Year), linetype = 2, color = "gray60") + # Optional: Group by year
      labs(title = paste(project_id, "- Community Shift by Year"),
           subtitle = "Shapes represent Locations, Colors represent Years",
           color = "Year",
           shape = "Location") +
      scale_color_brewer(palette = "Set1") # Gives distinct, high-contrast colors for years
    
    # 4. Save
    ggsave(file.path(figures_dir, paste0(project_id, "_pcoa_year_color_local_shape.png")), 
           p_beta_styled, width = 10, height = 7)
    
  }, error = function(e) {
    message("Warning: Could not generate the styled PCoA plot: ", e$message)
  })
  
  message("Generating PCoA (Jaccard): Color = Season, Shape = Location")
  
  tryCatch({
    # 1. Extract, Clean, and Re-insert Metadata
    # We extract to a standard data frame so dplyr works
    df_meta <- data.frame(sample_data(ps_object)) 
    
    df_meta <- df_meta %>%
      mutate(
        # Standardize the naming to match our color vector exactly
        season = case_when(
          grepl("non", season, ignore.case = TRUE) ~ "Non-upwelling",
          grepl("upwelling", season, ignore.case = TRUE) ~ "Upwelling",
          TRUE ~ as.character(season)
        ),
        # Use your helper for the locations
        local = clean_grouping(local)
      )
    
    # IMPORTANT: Put the rownames back (phyloseq needs these to match OTUs)
    rownames(df_meta) <- sample_names(ps_object)
    
    # Re-wrap as sample_data and update the phyloseq object
    sample_data(ps_object) <- sample_data(df_meta)
    
    # 2. Run Ordination
    ord_pcoa <- ordinate(ps_object, method = "PCoA", distance = "jaccard", binary = TRUE)
    
    # 3. Define Palette
    season_palette <- c(
      "Upwelling"     = "#56B4E9", 
      "Non-upwelling" = "#E69F00"
    )
    
    # 4. Plot
    p_season_local <- plot_ordination(ps_object, ord_pcoa, 
                                      color = "season", 
                                      shape = "local") +
      geom_point(size = 6, alpha = 0.8) +
      theme_minimal() +
      labs(title = paste(project_id, "- Seasonal & Spatial Community Structure"),
           color = "Season",
           shape = "Location") +
      scale_color_manual(values = season_palette)
    
    # Only add interaction ellipses if we have enough samples per combination
    interaction_counts <- table(interaction(df_meta$season, df_meta$local))
    if (any(interaction_counts >= 4)) {
      p_season_local <- p_season_local + 
        stat_ellipse(aes(group = interaction(season, local)), linetype = 2, na.rm = TRUE)
    } 
    
    # 5. Save
    ggsave(file.path(figures_dir, paste0(project_id, "_pcoa_season_location.png")), 
           p_season_local, width = 10, height = 7)
    
  }, error = function(e) {
    message("Warning: Could not generate Season/Location PCoA: ", e$message)
  })
  
}
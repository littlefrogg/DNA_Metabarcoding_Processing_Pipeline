generate_otu_barplots <- function(ps_object, figures_dir, project_id, tax_palette, color_config = NULL) {
  
  # --- 1. Load Packages ---
  suppressPackageStartupMessages({
    require(phyloseq)
    require(ggplot2)
    require(dplyr)
    require(forcats)
    require(tidyr)
  })
  
  # --- 2. Load Helper Functions ---
  source(file.path(dirname(figures_dir), "..", "scripts", "taxa_color_utils.R"))
  
  # --- 3. Data Cleaning & Metadata Standardization ---
  clean_grouping <- function(x) {
    x <- trimws(as.character(x))
    x <- tolower(x)
    x <- gsub("_", " ", x)
    x <- tools::toTitleCase(x)
    x
  }

  color_config <- resolve_color_config(color_config)
  
  # Extract and fix metadata
  df_meta <- data.frame(sample_data(ps_object)) %>%
    mutate(
      season = case_when(
        grepl("non", season, ignore.case = TRUE) ~ "Non-upwelling",
        grepl("upwelling", season, ignore.case = TRUE) ~ "Upwelling",
        TRUE ~ as.character(season)
      ),
      local = clean_grouping(local)
    )
  rownames(df_meta) <- sample_names(ps_object)
  sample_data(ps_object) <- sample_data(df_meta)
  
  if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
  
  # --- 4. Get Highest Taxonomic Resolution ---
  get_highest_taxonomy <- function(ps) {
    tax_mat <- as.data.frame(tax_table(ps))
    
    # For each OTU, find first non-empty value from right to left (most specific to least)
    highest_tax <- apply(tax_mat, 1, function(row) {
      for (i in length(row):1) {
        val <- row[i]
        if (!is.na(val) && val != "" && val != "unassigned" && tolower(val) != "no identification") {
          return(val)
        }
      }
      return("Unknown")
    })
    
    return(highest_tax)
  }
  
  # --- 5. Create OTU labels with highest taxonomy ---
  create_otu_labels <- function(ps) {
    highest_tax <- get_highest_taxonomy(ps)
    otu_names <- taxa_names(ps)
    
    # Create labels: OTU_name (highest_taxonomy)
    otu_labels <- paste0(otu_names, " (", highest_tax, ")")
    names(otu_labels) <- otu_names
    
    return(otu_labels)
  }
  
  # --- 5B. Create OTU labels with rank and highest taxonomy ---
  create_otu_labels_with_rank <- function(ps) {
    tax_mat <- as.data.frame(tax_table(ps))
    otu_names <- taxa_names(ps)
    
    # For each OTU, find first non-empty value from right to left and capture the rank name
    otu_labels <- apply(tax_mat, 1, function(row) {
      for (i in length(row):1) {
        val <- row[i]
        rank <- names(row)[i]
        if (!is.na(val) && val != "" && val != "unassigned" && tolower(val) != "no identification") {
          return(paste0(rank, ": ", val))
        }
      }
      return("Unassigned")
    })
    
    # Prepend OTU name to keep labels unique
    otu_labels <- paste0(otu_names, " (", otu_labels, ")")
    names(otu_labels) <- otu_names
    return(otu_labels)
  }
  
  # --- 6. Helper Function for Highest Rank Barplots ---
  save_otu_highest_rank_barplot <- function(ps, title_suffix, filename, x_var = "local", facet_var = NULL) {
    message("- Generating OTU highest rank barplot for: ", title_suffix)
    
    # Create a copy to avoid modifying original
    ps_copy <- ps
    
    # Get OTU labels with rank and highest taxonomy
    otu_labels_with_rank <- create_otu_labels_with_rank(ps_copy)
    
    # Rename taxa to include rank and highest taxonomy
    taxa_names(ps_copy) <- otu_labels_with_rank[taxa_names(ps_copy)]
    
    # Convert to presence/absence for richness counting
    otu_table(ps_copy) <- ifelse(otu_table(ps_copy) > 0, 1, 0)
    
    # Add columns to tax_table with the highest rank value and rank name
    tax_mat <- as.data.frame(tax_table(ps))
    
    # Extract the highest rank value and rank name for each OTU
    highest_rank_values <- apply(tax_mat, 1, function(row) {
      for (i in length(row):1) {
        val <- row[i]
        if (!is.na(val) && val != "" && val != "unassigned" && tolower(val) != "no identification") {
          return(val)
        }
      }
      return("Unassigned")
    })
    
    highest_rank_names <- apply(tax_mat, 1, function(row) {
      for (i in length(row):1) {
        val <- row[i]
        rank <- names(row)[i]
        if (!is.na(val) && val != "" && val != "unassigned" && tolower(val) != "no identification") {
          return(rank)
        }
      }
      return("Unassigned")
    })
    
    # Add as new columns
    tax_mat_copy <- as.data.frame(tax_table(ps_copy))
    tax_mat_copy$HighestTaxon <- highest_rank_values
    tax_mat_copy$HighestRank <- highest_rank_names
    
    # Convert HighestTaxon to a factor ordered by Phylum and HighestTaxon for proper stacking
    ordered_factor_levels <- tax_mat_copy %>%
      arrange(Phylum, HighestTaxon) %>%
      pull(HighestTaxon) %>%
      unique()
    
    tax_mat_copy$HighestTaxon <- factor(tax_mat_copy$HighestTaxon, levels = ordered_factor_levels)
    
    tax_table(ps_copy) <- as.matrix(tax_mat_copy)
    
    # Build color mapping: use base palette for Phylum, use shades of phylum colors for lower ranks
    # Map taxa to colors based on rank and phylum
    my_colors_map <- list()
    
    # For Phylum: use base palette
    if ("Phylum" %in% highest_rank_names) {
      phylum_taxa <- unique(highest_rank_values[highest_rank_names == "Phylum"])
      for (taxon in phylum_taxa) {
        if (taxon %in% names(tax_palette)) {
          my_colors_map[[taxon]] <- tax_palette[taxon]
        }
      }
    }
    
    # For non-Phylum ranks: generate shades of their phylum's base color
    for (rank in c("Class", "Order", "Family", "Genus", "Species")) {
      if (rank %in% highest_rank_names) {
        # Get all OTUs at this rank and their phyla
        rank_indices <- which(highest_rank_names == rank)
        
        # Group by phylum within this rank
        rank_by_phylum <- list()
        for (idx in rank_indices) {
          phylum <- tax_mat_copy[idx, "Phylum"]
          taxon_val <- highest_rank_values[idx]
          
          if (!(taxon_val %in% names(my_colors_map))) {
            if (!(phylum %in% names(rank_by_phylum))) {
              rank_by_phylum[[phylum]] <- c()
            }
            rank_by_phylum[[phylum]] <- c(rank_by_phylum[[phylum]], taxon_val)
          }
        }
        
        # For each phylum group at this rank, generate shades
        for (phylum in names(rank_by_phylum)) {
          phylum_taxa_at_rank <- unique(rank_by_phylum[[phylum]])
          
          if (phylum %in% names(tax_palette)) {
            base_col <- tax_palette[phylum]
            # Generate shades from lighter to darker through base
            ramp <- colorRampPalette(
              c(colorspace::lighten(base_col, color_config$highest_rank_light_low),
                colorspace::lighten(base_col, color_config$highest_rank_light_high))
            )(length(phylum_taxa_at_rank))
            
            for (i in seq_along(phylum_taxa_at_rank)) {
              my_colors_map[[phylum_taxa_at_rank[i]]] <- ramp[i]
            }
          }
        }
      }
    }
    
    # Convert to named vector
    my_colors_list <- sapply(highest_rank_values, function(taxon) {
      if (taxon %in% names(my_colors_map)) {
        return(my_colors_map[[taxon]])
      } else {
        return("#CCCCCC")
      }
    }, simplify = TRUE)
    
    # Convert to named vector
    my_colors <- setNames(my_colors_list, highest_rank_values)
    my_colors <- my_colors[!duplicated(names(my_colors))]
    
    # Reorder legend by Phylum for a professional look
    ordered_levels <- tax_mat_copy %>%
      arrange(Phylum, HighestTaxon) %>%
      pull(HighestTaxon) %>%
      unique()
    
    # Create Plot using manual ggplot for proper stacking control
    plot_data <- psmelt(ps_copy)
    plot_data$HighestTaxon <- factor(plot_data$HighestTaxon, levels = ordered_levels)
    
    p <- ggplot(plot_data, aes(x = !!sym(x_var), y = Abundance, fill = HighestTaxon)) +
      geom_bar(stat = "identity", position = "stack", alpha = color_config$bar_alpha_highest_rank, color = NA) +
      scale_fill_manual(values = my_colors, breaks = ordered_levels) +
      theme_minimal() +
      labs(title = paste(project_id, "-", title_suffix),
           y = "Richness (Number of Unique OTUs)", 
           x = tools::toTitleCase(x_var),
           fill = "Highest Taxon") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.key.size = unit(0.3, "cm"),
            legend.text = element_text(size = 6),
            strip.background = element_rect(fill = "gray90"),
            strip.text = element_text(face = "bold"))
    
    # Add faceting if specified
    if (!is.null(facet_var)) {
      p <- p + facet_wrap(as.formula(paste("~", facet_var)), scales = "free_x")
    }
    
    ggsave(file.path(figures_dir, filename), p, width = 14, height = 10)
  }
  
  # --- 7. Helper Function to Create regular OTU Barplots with proper stacking ---
  save_otu_barplot_with_stacking <- function(ps, title_suffix, filename, tax_level = "Phylum", x_var = "local", facet_var = NULL) {
    message("- Generating OTU barplot with stacking for: ", title_suffix)
    
    # Get OTU labels with highest taxonomy
    otu_labels <- create_otu_labels(ps)
    
    # Rename taxa to include highest taxonomy
    taxa_names(ps) <- otu_labels[taxa_names(ps)]
    
    # Convert to presence/absence for richness counting
    otu_table(ps) <- ifelse(otu_table(ps) > 0, 1, 0)
    
    # Get consistent colors using centralized utilities
    is_phylum_plot <- tolower(tax_level) == "phylum"
    my_colors <- get_plot_taxa_colors(ps, level = tax_level, base_palette = tax_palette, color_config = color_config)
    
    # Reorder legend by Phylum for a professional look
    tax_df <- as.data.frame(tax_table(ps))
    ordered_levels <- tax_df %>%
      arrange(Phylum, !!sym(tax_level)) %>%
      pull(!!sym(tax_level)) %>%
      unique()
    
    # Create Plot using manual ggplot for proper stacking control
    plot_data <- psmelt(ps)
    plot_data[[tax_level]] <- factor(plot_data[[tax_level]], levels = ordered_levels)
    
    p <- ggplot(plot_data, aes(x = !!sym(x_var), y = Abundance, fill = !!sym(tax_level))) +
      geom_bar(stat = "identity", position = "stack", alpha = if (is_phylum_plot) color_config$bar_alpha_phylum else color_config$bar_alpha_other, color = NA) +
      scale_fill_manual(values = my_colors, breaks = ordered_levels) +
      theme_minimal() +
      labs(title = paste(project_id, "-", title_suffix),
           y = "Richness (Number of Unique OTUs)", 
           x = tools::toTitleCase(x_var)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.key.size = unit(0.3, "cm"),
            legend.text = element_text(size = 7),
            strip.background = element_rect(fill = "gray90"),
            strip.text = element_text(face = "bold"))
    
    # Add faceting if specified
    if (!is.null(facet_var)) {
      p <- p + facet_wrap(as.formula(paste("~", facet_var)), scales = "free_x")
    }
    
    ggsave(file.path(figures_dir, filename), p, width = 14, height = 10)
  }
  # --- 8. Shared Helper: run plot suite for one taxonomic subset ---
  run_season_local_plot_suite <- function(ps_subset, prefix, label) {
    sample_data(ps_subset)$MergeKey <- paste(sample_data(ps_subset)$season,
                                             sample_data(ps_subset)$local, sep = "___")
    ps_merged <- merge_samples(ps_subset, "MergeKey")

    new_meta <- data.frame(MergeKey = rownames(sample_data(ps_merged))) %>%
      separate(MergeKey, into = c("season", "local"), sep = "___", remove = FALSE)
    rownames(new_meta) <- new_meta$MergeKey
    sample_data(ps_merged) <- sample_data(new_meta)

    save_otu_barplot_with_stacking(ps_merged, paste(label, "OTUs by Season & Location"),
                                   paste0(project_id, "_", prefix, "_otus_seasonal_faceted.png"),
                                   tax_level = "Phylum", x_var = "season", facet_var = "local")

    save_otu_barplot_with_stacking(ps_merged, paste(label, "OTUs (Class) by Season & Location"),
                                   paste0(project_id, "_", prefix, "_otus_class_seasonal_faceted.png"),
                                   tax_level = "Class", x_var = "season", facet_var = "local")

    save_otu_barplot_with_stacking(ps_merged, paste(label, "OTUs (Order) by Season & Location"),
                                   paste0(project_id, "_", prefix, "_otus_order_seasonal_faceted.png"),
                                   tax_level = "Order", x_var = "season", facet_var = "local")

    save_otu_barplot_with_stacking(ps_merged, paste(label, "OTUs (Family) by Season & Location"),
                                   paste0(project_id, "_", prefix, "_otus_family_seasonal_faceted.png"),
                                   tax_level = "Family", x_var = "season", facet_var = "local")

    save_otu_barplot_with_stacking(ps_merged, paste(label, "OTUs (Genus) by Season & Location"),
                                   paste0(project_id, "_", prefix, "_otus_genus_seasonal_faceted.png"),
                                   tax_level = "Genus", x_var = "season", facet_var = "local")

    save_otu_barplot_with_stacking(ps_merged, paste(label, "OTUs (Species) by Season & Location"),
                                   paste0(project_id, "_", prefix, "_otus_species_seasonal_faceted.png"),
                                   tax_level = "Species", x_var = "season", facet_var = "local")

    save_otu_highest_rank_barplot(ps_merged, paste(label, "OTUs (Highest Rank) by Season & Location"),
                                  paste0(project_id, "_", prefix, "_otus_highest_rank_seasonal_faceted.png"),
                                  x_var = "season", facet_var = "local")
  }

  target_phyla <- c("Arthropoda", "Mollusca", "Annelida", "Nemertea", 
                    "Cnidaria", "Porifera", "Echinodermata", "Gastrotricha")
  tt <- as.data.frame(tax_table(ps_object))
  pattern <- paste0("^(", paste(target_phyla, collapse = "|"), ")(_|$)")
  keep_taxa <- rownames(tt)[grepl(pattern, tt$Phylum, ignore.case = TRUE)]
  
  if (length(keep_taxa) > 0) {
    ps_invert <- prune_taxa(keep_taxa, ps_object)
    run_season_local_plot_suite(ps_invert, prefix = "invertebrates", label = "Invertebrate")
  }
  
  # --- 7. Process Chordata ---
  chordata_pattern <- "^(Chordata)(_|$)"
  keep_chordata <- rownames(tt)[grepl(chordata_pattern, tt$Phylum, ignore.case = TRUE)]
  
  if (length(keep_chordata) > 0) {
    ps_chordata <- prune_taxa(keep_chordata, ps_object)
    run_season_local_plot_suite(ps_chordata, prefix = "chordata", label = "Chordata")
  } else {
    message("Warning: No Chordata taxa found.")
  }
  
  message("OTU barplots generated successfully.")
}
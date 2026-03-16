generate_taxa_plots <- function(ps_object, figures_dir, project_id, tax_palette) {
  
  # --- 1. Load Packages ---
  suppressPackageStartupMessages({
    require(phyloseq)
    require(ggplot2)
    require(dplyr)
    require(forcats)
    require(tidyr)
  })
  
  # --- 2. Data Cleaning & Metadata Standardization ---
  clean_grouping <- function(x) {
    x <- trimws(as.character(x))
    x <- tolower(x)
    x <- gsub("_", " ", x)
    x <- tools::toTitleCase(x)
    x
  }
  
  # Extract and fix metadata to ensure Season colors match the palette
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
  
  # --- 3. Robust Plotting Helper ---
  # Handles x-axis selection (Local or Season) and optional faceting
  save_taxa_plot <- function(ps, level, title_suffix, filename, x_var = "local", facet_var = NULL) {
    message("- Generating richness plot for: ", title_suffix)
    
    # Agglomerate to target level
    ps_glom <- tax_glom(ps, level)
    
    # Get consistent colors from your utility script
    my_colors <- get_taxa_colors(ps_glom, level = level, base_palette = tax_palette)
    
    # Reorder legend by Phylum for a professional look
    tax_df <- as.data.frame(tax_table(ps_glom))
    ordered_levels <- tax_df %>%
      arrange(Phylum, !!sym(level)) %>%
      pull(!!sym(level))
    
    # Convert to presence/absence for richness (Count = 1 per OTU)
    otu_table(ps_glom) <- ifelse(otu_table(ps_glom) > 0, 1, 0)
    
    # Create Plot
    p <- plot_bar(ps_glom, x = x_var, fill = level) +
      geom_bar(aes(color = !!sym(level), fill = !!sym(level)), 
               stat = "identity", position = "stack") +
      scale_fill_manual(values = my_colors, breaks = ordered_levels) +
      scale_color_manual(values = my_colors, breaks = ordered_levels) +
      theme_minimal() +
      labs(title = paste(project_id, "-", title_suffix),
           y = "Richness (Number of OTUs)", 
           x = tools::toTitleCase(x_var)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.key.size = unit(0.3, "cm"),
            legend.text = element_text(size = 8),
            strip.background = element_rect(fill = "gray90"), # Facet label background
            strip.text = element_text(face = "bold"))
    
    # Add faceting if a variable (like 'local') is provided
    if (!is.null(facet_var)) {
      p <- p + facet_wrap(as.formula(paste("~", facet_var)), scales = "free_x")
    }
    
    ggsave(file.path(figures_dir, filename), p, width = 14, height = 8)
  }
  
  # --- 4. Process Invertebrates ---
  target_phyla <- c("Arthropoda", "Mollusca", "Annelida", "Nemertea", 
                    "Cnidaria", "Porifera", "Echinodermata", "Gastrotricha")
  tt <- as.data.frame(tax_table(ps_object))
  pattern <- paste0("^(", paste(target_phyla, collapse = "|"), ")(_|$)")
  keep_taxa <- rownames(tt)[grepl(pattern, tt$Phylum, ignore.case = TRUE)]
  
  if (length(keep_taxa) > 0) {
    ps_invert <- prune_taxa(keep_taxa, ps_object)
    
    # --- Prepare Data for Faceted Plots ---
    # merge_samples clears out metadata, so we create a combined key to recover it
    sample_data(ps_invert)$MergeKey <- paste(sample_data(ps_invert)$season, 
                                             sample_data(ps_invert)$local, sep = "___")
    
    ps_merged <- merge_samples(ps_invert, "MergeKey")
    
    # Repair metadata columns after merge
    new_meta <- data.frame(MergeKey = rownames(sample_data(ps_merged))) %>%
      separate(MergeKey, into = c("season", "local"), sep = "___", remove = FALSE)
    rownames(new_meta) <- new_meta$MergeKey
    sample_data(ps_merged) <- sample_data(new_meta)

# A. Invertebrate Richness by Phylum (Season, faceted by Location)
    save_taxa_plot(ps_merged, "Phylum", "Invertebrate Richness (Phylum) by Season & Location", 
                   paste0(project_id, "_invertebrates_phylum_seasonal_faceted.png"), 
                   x_var = "season", facet_var = "local")
    
# B. Invertebrate Richness by Class (Season, faceted by Location)
    save_taxa_plot(ps_merged, "Class", "Invertebrate Richness (Class) by Season & Location", 
                   paste0(project_id, "_invertebrates_class_seasonal_faceted.png"), 
                   x_var = "season", facet_var = "local")
    
# C. Invertebrate Richness by Season (Faceted by Location)
    save_taxa_plot(ps_merged, "Order", "Invertebrate Richness (Order) by Season & Location", 
                   paste0(project_id, "_invertebrates_order_seasonal_faceted.png"), 
                   x_var = "season", facet_var = "local")
    
# D. Invertebrate Genus Level by Season (Faceted by Location)
    save_taxa_plot(ps_merged, "Genus", "Invertebrate Richness (Genus) by Season", 
                   paste0(project_id, "_invertebrates_genus_seasonal_faceted.png"), 
                   x_var = "season", facet_var = "local")
    
# E. Invertebrate Species Level by Season (Faceted by Location)
    save_taxa_plot(ps_merged, "Species", "Invertebrate Richness (Species) by Season", 
                   paste0(project_id, "_invertebrates_species_seasonal_faceted.png"), 
                   x_var = "season", facet_var = "local")

  } else {
    message("Warning: No taxa found for the specified invertebrate phyla.")
  }
  
  # --- 5. Process Chordata ---
  chordata_pattern <- "^(Chordata)(_|$)"
  keep_chordata <- rownames(tt)[grepl(chordata_pattern, tt$Phylum, ignore.case = TRUE)]
  
  if (length(keep_chordata) > 0) {
    ps_chordata <- prune_taxa(keep_chordata, ps_object)
    
    # Same merge logic for Chordata
    sample_data(ps_chordata)$MergeKey <- paste(sample_data(ps_chordata)$season, 
                                               sample_data(ps_chordata)$local, sep = "___")
    ps_chord_merged <- merge_samples(ps_chordata, "MergeKey")
    
    new_meta_c <- data.frame(MergeKey = rownames(sample_data(ps_chord_merged))) %>%
      separate(MergeKey, into = c("season", "local"), sep = "___", remove = FALSE)
    rownames(new_meta_c) <- new_meta_c$MergeKey
    sample_data(ps_chord_merged) <- sample_data(new_meta_c)
    
    save_taxa_plot(ps_chord_merged, "Order", "Chordata Richness (Order) by Season & Location", 
                   paste0(project_id, "_chordata_order_seasonal_faceted.png"), 
                   x_var = "season", facet_var = "local")
  } else {
    message("Warning: No Chordata taxa found.")
  }
  
  message("Additional plots generated successfully.")
}
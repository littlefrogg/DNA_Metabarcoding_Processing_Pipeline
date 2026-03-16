################################################################
# GET TAXA COLORS
################################################################
# Purpose: This script defines a function to assign colors to OTUS
#
# Author: Paige Smallman, 2026
################################################################
# scripts/taxa_color_utils.R

#' Generate Consistent Taxonomic Colors
#' @param ps A phyloseq object
#' @param level The taxonomic level to color (e.g., "Phylum", "Class", "Genus")
#' @param base_palette A named vector of colors (e.g., Phylum_palette)
#' @param default_color Hex code for phyla not in the palette
#' @param na_color Hex code for "no identification"
get_taxa_colors <- function(ps, level = "Phylum", base_palette, 
                            default_color = "#1300a0ff", 
                            na_color = "#d0d0d0ff") {
  require(dplyr)
  require(tibble)
  require(colorspace)
  
  # 1. Extract tax table and clean suffixes (Arthropoda_6656 -> Arthropoda)
  tax_df <- as.data.frame(tax_table(ps)@.Data) %>%
    mutate(across(everything(), ~ gsub("_\\d+", "", .x)))
  
  # 2. Map unique taxa at the target level to their Phylum
  tax_map <- tax_df %>%
    select(Phylum, !!sym(level)) %>%
    distinct()
  
  unique_taxa <- unique(tax_map[[level]])
  
  # Initialize with default color
  final_colors <- setNames(rep(default_color, length(unique_taxa)), unique_taxa)
  
  # 3. Generate shades for each Phylum in the base_palette
  for (phy in names(base_palette)) {
    sub_taxa <- tax_map %>% 
      filter(Phylum == phy) %>% 
      pull(!!sym(level)) %>% 
      unique()
    
    # Exclude "no identification" from the color ramp
    sub_taxa_id <- sub_taxa[sub_taxa != "no identification"]
    
    if (length(sub_taxa_id) > 0) {
      base_col <- base_palette[phy]
      # Create shades from darker to lighter through base color
      ramp <- colorRampPalette(c(darken(base_col, 0.5), base_col, lighten(base_col, 0.7)))(length(sub_taxa_id))
      final_colors[sub_taxa_id] <- ramp
    }
  }
  
  # 4. Final override for "no identification"
  if ("no identification" %in% names(final_colors)) {
    final_colors["no identification"] <- na_color
  }
  
  return(final_colors)
}

#' Default plotting color configuration
get_default_color_config <- function() {
  list(
    use_base_for_phylum = TRUE,
    lighten_non_phylum = 0.18,
    bar_alpha_phylum = 1.00,
    bar_alpha_other = 0.85,
    bar_alpha_highest_rank = 0.88,
    highest_rank_light_low = 0.45,
    highest_rank_light_high = 0.15,
    fallback_color = "#BDBDBD"
  )
}

#' Merge user config over defaults
resolve_color_config <- function(color_config = NULL) {
  defaults <- get_default_color_config()
  if (is.null(color_config)) return(defaults)
  defaults[names(color_config)] <- color_config
  defaults
}

#' Get exact base phylum colors (no shading)
get_base_phylum_colors <- function(ps, level = "Phylum", base_palette, fallback_color = "#BDBDBD") {
  tax_df <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
  taxa_vals <- as.character(tax_df[[level]])
  taxa_vals[is.na(taxa_vals) | trimws(taxa_vals) == ""] <- "Unknown"
  taxa_clean <- gsub("_[0-9]+$", "", taxa_vals)
  cols <- ifelse(taxa_clean %in% names(base_palette), base_palette[taxa_clean], fallback_color)
  cols <- setNames(as.character(cols), taxa_vals)
  cols[!duplicated(names(cols))]
}

#' Lighten a named color vector by amount
lighten_color_vector <- function(cols, amount = 0.18) {
  vals <- vapply(cols, function(cl) colorspace::lighten(cl, amount), character(1))
  setNames(vals, names(cols))
}

#' Central color getter for OTU/taxa barplots
get_plot_taxa_colors <- function(ps, level = "Phylum", base_palette, color_config = NULL) {
  cfg <- resolve_color_config(color_config)
  is_phylum <- tolower(level) == "phylum"

  if (is_phylum && isTRUE(cfg$use_base_for_phylum)) {
    return(get_base_phylum_colors(ps, level = level, base_palette = base_palette, fallback_color = cfg$fallback_color))
  }

  cols <- get_taxa_colors(ps, level = level, base_palette = base_palette)
  if (!is.null(cfg$lighten_non_phylum) && cfg$lighten_non_phylum > 0) {
    cols <- lighten_color_vector(cols, amount = cfg$lighten_non_phylum)
  }
  cols
}
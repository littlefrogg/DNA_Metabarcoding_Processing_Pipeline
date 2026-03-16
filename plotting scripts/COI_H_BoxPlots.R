################################################################
# (COI-H) BOX PLOTS - DNA metabarcoding detailed visualization
################################################################
# Purpose: Generate additional box plots for the final phyloseq object.
#
# Author: Paige Smallman, 2026
################################################################
generate_alpha_boxplots <- function(ps_object, figures_dir, project_id, permanova_p = NULL, permanova_r2 = NULL, pairwise_results = NULL) {
	# Define season palette
	       season_palette <- c(
		       "Upwelling"     = "#56B4E9",
		       "Non-Upwelling" = "#E69F00"
	       )
suppressPackageStartupMessages({
	require(phyloseq)
	require(ggplot2)
	require(dplyr)
	require(tidyr)
	if (!is.null(pairwise_results)) require(ggsignif)
	})

	# Clean grouping columns
	clean_grouping <- function(x) {
		x <- trimws(as.character(x))
		x <- tolower(x)
		x <- gsub("_", " ", x)
		x <- tools::toTitleCase(x)
		x
	}
	       df_meta <- data.frame(sample_data(ps_object)) %>%
		       mutate(
			       season = clean_grouping(season),
			       local = clean_grouping(local),
			       year = if("year" %in% names(.)) clean_grouping(year) else if("Year" %in% names(.)) clean_grouping(Year) else NA
		       )
	df_meta$Sample <- rownames(df_meta)
	rownames(df_meta) <- sample_names(ps_object)
	sample_data(ps_object) <- sample_data(df_meta)

	if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

# Calculate alpha diversity (Observed)
	alpha_df <- estimate_richness(ps_object, measures = "Observed")
	alpha_df <- alpha_df %>%
		rownames_to_column("Sample") %>%
		left_join(df_meta, by = "Sample")

		# Add subtitle if permanova_p and permanova_r2 are provided
		subtitle_text <- NULL
	#	if (!is.null(permanova_p) && !is.null(permanova_r2) && !is.na(permanova_p) && !is.na(permanova_r2)) {
	#		subtitle_text <- paste0("Upwelling vs. Non-Upwelling (both locations) PERMANOVA: p = ", formatC(as.numeric(permanova_p), format = "e", digits = 2), ", R² = ", round(as.numeric(permanova_r2), 3))
	#	}
		p1 <- ggplot(alpha_df, aes(x = season, y = Observed, fill = season)) +
				geom_boxplot(alpha = 0.7) +
				facet_wrap(~ local, scales = "free_x") +
				scale_fill_manual(values = season_palette) +
				theme_bw() +
				labs(title = paste(project_id, "- Invertebrate Alpha Diversity by Season and Location"),
						 subtitle = subtitle_text,
						 x = "Season", y = "Observed Richness") +
				theme(axis.text.x = element_text(angle = 45, hjust = 1))

	ggsave(file.path(figures_dir, paste0(project_id, "_alpha_boxplot_season_faceted_location.png")), p1, width = 12, height = 7)
	message("Alpha diversity boxplot (season, faceted by location) saved.")


# Create a year_season variable for x-axis and facet by local
			 year_season_levels <- c(
				 "2020_Non-Upwelling",
				 "2021_Upwelling",
				 "2021_Non-Upwelling",
				 "2022_Upwelling",
				 "2022_Non-Upwelling",
				 "2023_Upwelling",
				 "2023_Non-Upwelling",
				 "2024_Upwelling",
				 "2024_Non-Upwelling",
				 "2025_Upwelling"
			 )
			 sample_data(ps_object)$year_season <- factor(
				 paste(sample_data(ps_object)$Year, clean_grouping(sample_data(ps_object)$season), sep = "_"),
				 levels = year_season_levels
			 )
			 sample_data(ps_object)$local <- clean_grouping(sample_data(ps_object)$local)
	
	# Add ENSO column based on date ranges
	if ("date_col" %in% colnames(df_meta)) {
	  
	  # --- TROUBLESHOOTING: Inspect date conversion ---
	  message("\n--- DEBUGGING ENSO PLOT ---")
	  message("1. Original 'date_col' (first 6 rows):")
	  message(paste(head(df_meta$date_col), collapse = ", "))
	  
	  df_meta <- df_meta %>%
	    mutate(
	      Date = as.Date(date_col, format = "%d-%m-%Y"), # Corrected format
	      ENSO = case_when(
	        Date >= as.Date("2020-01-01") & Date <= as.Date("2020-08-31") ~ "Neutral (2020)",
	        Date >= as.Date("2020-08-01") & Date <= as.Date("2021-05-31") ~ "La Niña (2020-21)",
	        Date >= as.Date("2021-05-01") & Date <= as.Date("2021-08-31") ~ "Neutral (2021)",
	        Date >= as.Date("2021-08-01") & Date <= as.Date("2023-01-31") ~ "La Niña (2021-23)",
	        Date >= as.Date("2023-01-01") & Date <= as.Date("2023-05-31") ~ "Neutral (2023)",
	        Date >= as.Date("2023-05-01") & Date <= as.Date("2024-04-30") ~ "El Niño (2023-24)",
	        Date >= as.Date("2024-05-01") & Date <= as.Date("2026-01-31") ~ "Neutral (2024-26)",
	        TRUE ~ "Unknown"
	      )
	    )
	    
	  message("\n2. Converted 'Date' column (first 6 rows):")
	  message(paste(head(df_meta$Date), collapse = ", "))
	  message("\n3. Resulting 'ENSO' column (first 6 rows):")
	  message(paste(head(df_meta$ENSO), collapse = ", "))
	  
	  # --- NEW: Inspect rows that are "Unknown" ---
	  unknown_enso_dates <- df_meta %>%
	    filter(ENSO == "Unknown") %>%
	    select(Sample, date_col, Date)
	  
	  if (nrow(unknown_enso_dates) > 0) {
	    message("\n4. Dates that were not assigned an ENSO phase:")
	    print(unknown_enso_dates)
	  } else {
	    message("\n4. All dates were successfully assigned an ENSO phase.")
	  }
	  
	  message("--- END DEBUGGING ---\n")
	  
	  # Re-join with alpha_df
	  alpha_df <- alpha_df %>%
	    select(-any_of(c("ENSO", "Date.y"))) %>% # Remove old ENSO if it exists
	    left_join(select(df_meta, Sample, ENSO, Date), by = "Sample")
	  
	  # Create ENSO plot
	  
	  # Define ENSO color palettes
	  enso_palette <- setNames(
	    unlist(lapply(unique(alpha_df$ENSO), function(phase) {
	      if (grepl("El Niño", phase)) return("#dc3620ff")
	      if (grepl("La Niña", phase)) return("#057abdff")
	      if (grepl("Neutral", phase)) return("#999999")
	      return("black") # Fallback for "Unknown"
	    })),
	    unique(alpha_df$ENSO)
	  )
	  
	  p_enso <- ggplot(alpha_df, aes(x = ENSO, y = Observed, fill = ENSO)) +
	    geom_boxplot(alpha = 0.7) +
	    facet_wrap(~local, scales = "free_x") +
	    scale_fill_manual(values = enso_palette) +
	    theme_bw() +
	    labs(title = paste(project_id, "- Invertebrate Alpha Diversity by ENSO Phase and Location"),
	         x = "ENSO Phase", y = "Observed Richness") +
	    theme(axis.text.x = element_text(angle = 45, hjust = 1))
	  
	  ggsave(file.path(figures_dir, paste0(project_id, "_alpha_enso_faceted.png")), p_enso, width = 16, height = 7)
	  message("Alpha diversity boxplot (ENSO, faceted by location) saved.")
	  
	  # --- NEW: Create aggregated ENSO plot ---
	  # Create a simplified ENSO phase column
	  alpha_df <- alpha_df %>%
	    mutate(ENSO_simple = gsub(" \\(.*\\)", "", ENSO))
	    
	  # Define order for the simplified phases
	  enso_simple_levels <- c("El Niño", "La Niña", "Neutral")
	  alpha_df$ENSO_simple <- factor(alpha_df$ENSO_simple, levels = intersect(enso_simple_levels, unique(alpha_df$ENSO_simple)))
	  
	  # Define simple palette
	  enso_simple_palette <- c(
	    "El Niño" = "#dc3620ff",
	    "La Niña" = "#057abdff",
	    "Neutral" = "#999999"
	  )
	  
	  # Create the aggregated plot
	  p_enso_simple <- ggplot(alpha_df, aes(x = ENSO_simple, y = Observed, fill = ENSO_simple)) +
	    geom_boxplot(alpha = 0.7) +
	    facet_wrap(~local, scales = "free_x") +
	    scale_fill_manual(values = enso_simple_palette) +
	    theme_bw() +
	    labs(title = paste(project_id, "- Invertebrate Alpha Diversity by Aggregated ENSO Phase"),
	         x = "ENSO Phase", y = "Observed Richness") +
	    theme(axis.text.x = element_text(angle = 45, hjust = 1))
	    
	  ggsave(file.path(figures_dir, paste0(project_id, "_alpha_enso_simple_faceted.png")), p_enso_simple, width = 12, height = 7)
	  message("Alpha diversity boxplot (Aggregated ENSO, faceted by location) saved.")
	  
	} else {
	  message("Skipping ENSO plot: 'date_col' column not found in metadata.")
	}
	
p2 <- plot_richness(ps_object, x = "year_season", measures = c("Observed")) +
				geom_boxplot(aes(fill = season), alpha = 0.7) +
				scale_fill_manual(values = season_palette, drop = FALSE) +
				facet_wrap(~local, scales = "free_x") +
				theme_bw() +
				labs(title = paste(project_id, "- Invertebrate Alpha Diversity by Year, Season, and Location"),
					x = "Year, Season", y = "Observed Richness") +
				scale_x_discrete(labels = function(x) gsub("_", ", ", x)) +
				theme(axis.text.x = element_text(angle = 45, hjust = 1))
	ggsave(file.path(figures_dir, paste0(project_id, "_alpha_year_season_local_faceted.png")), p2, width = 16, height = 7)
	message("Alpha diversity boxplot (year & season, faceted by location) saved.")
	

# Boxplot grouped by year, colored by season, faceted by location
		p3 <- ggplot(alpha_df, aes(x = Year, y = Observed, fill = year)) +
			geom_boxplot(alpha = 0.7) +
			facet_wrap(~ local, scales = "free_x") +
			theme_bw() +
			labs(title = paste(project_id, "- Invertebrate Alpha Diversity by Year and Season (Faceted by Location)"),
					 x = "Year", y = "Observed Richness", fill = "Season") +
			theme(axis.text.x = element_text(angle = 45, hjust = 1))
	
	# Add significance brackets from pairwise PERMANOVA if available
	if (!is.null(pairwise_results) && length(pairwise_results) > 0) {
	  # DEBUG: Print information about pairwise results
	  message("\n=== DEBUGGING PAIRWISE RESULTS ===")
	  message("Number of locations in pairwise_results: ", length(pairwise_results))
	  message("Location names in pairwise_results: ", paste(names(pairwise_results), collapse = ", "))
	  message("Unique locations in alpha_df: ", paste(unique(alpha_df$local), collapse = ", "))
	  message("Unique years in alpha_df: ", paste(unique(alpha_df$Year), collapse = ", "))
	  
	  sig_annotations <- list()
	  y_max <- max(alpha_df$Observed, na.rm = TRUE)
	  bracket_height <- y_max * 0.05
	  current_y <- y_max + bracket_height
	  sig_count <- 0
	  
	  for (loc_name in names(pairwise_results)) {
	    pair_list <- pairwise_results[[loc_name]]
	    # pair_list is a pwadstrata object (which is a list of comparisons)
	    # Each element is a data frame with comparison results
	    
	    if (is.list(pair_list)) {
	      # Iterate through each comparison in the list
	      for (comp_name in names(pair_list)) {
	        comp_result <- pair_list[[comp_name]]
	        
	        # Check if this is a data frame (comparison result)
	        if (is.data.frame(comp_result) && "Pr(>F)" %in% colnames(comp_result)) {
	          # Get p-value from Model row
	          if ("Model" %in% rownames(comp_result)) {
	            p_val <- comp_result["Model", "Pr(>F)"]
	            
	            # If significant, extract the year pair
	            if (!is.na(p_val) && p_val < 0.05) {
	              # comp_name format is like "2024_vs_2021"
	              years <- as.character(unlist(strsplit(comp_name, "_vs_")))
	              if (length(years) == 2) {
	                message("  Found significant pair: ", comp_name, " (p=", round(p_val, 4), ") in ", loc_name)
	                annotation_key <- paste(loc_name, comp_name, sep = ":")
	                sig_annotations[[annotation_key]] <- list(
	                  location = loc_name,
	                  year1 = years[1],
	                  year2 = years[2],
	                  p_value = p_val,
	                  y_pos = current_y
	                )
	                current_y <- current_y + bracket_height
	                sig_count <- sig_count + 1
	              }
	            }
	          }
	        }
	      }
	    }
	  }
	  
	  # Add geom_signif layers for each significant pair by location
	  message("Total significant pairs found: ", sig_count)
	  message("Total annotations to add: ", length(sig_annotations))
	  
	  if (length(sig_annotations) > 0) {
	    message("Adding significance brackets to plot...")
	    
	    # Convert annotations to a data frame for easier manipulation
	    bracket_data <- data.frame()
	    for (annotation in sig_annotations) {
	      message("  Adding bracket: ", annotation$year1, " vs ", annotation$year2, " in ", annotation$location)
	      
	      # Get numeric x positions for the years
	      year_levels <- unique(sort(as.numeric(alpha_df$Year)))
	      x1 <- which(year_levels == as.numeric(annotation$year1))
	      x2 <- which(year_levels == as.numeric(annotation$year2))
	      
	      bracket_data <- rbind(bracket_data, data.frame(
	        local = annotation$location,
	        x1 = x1,
	        x2 = x2,
	        y = annotation$y_pos,
	        p_val = annotation$p_value,
	        stringsAsFactors = FALSE
	      ))
	    }
	    
	    # Add bracket lines and significance markers
	    p3 <- p3 + 
	      # Horizontal lines at top of bracket
	      geom_segment(
	        data = bracket_data,
	        aes(x = x1, xend = x2, y = y, yend = y),
	        inherit.aes = FALSE,
	        color = "black",
	        linewidth = 0.5
	      ) +
	      # Vertical lines at ends of bracket
	      geom_segment(
	        data = bracket_data,
	        aes(x = x1, xend = x1, y = y * 0.99, yend = y),
	        inherit.aes = FALSE,
	        color = "black",
	        linewidth = 0.5
	      ) +
	      geom_segment(
	        data = bracket_data,
	        aes(x = x2, xend = x2, y = y * 0.99, yend = y),
	        inherit.aes = FALSE,
	        color = "black",
	        linewidth = 0.5
	      ) +
	      # Significance stars
	      geom_text(
	        data = bracket_data,
	        aes(x = (x1 + x2) / 2, y = y * 1.02, label = ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", "*"))),
	        inherit.aes = FALSE,
	        size = 4,
	        vjust = -0.5
	      )
	  } else {
	    message("No significant pairs found to plot.")
	  }
	}
}
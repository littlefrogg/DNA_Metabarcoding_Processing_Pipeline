################################################################
# (COI-H) PERMANOVA ANALYSIS - DNA metabarcoding
################################################################
# Purpose: Run PERMANOVA to compare Upwelling vs Non-Upwelling (overall and within locations)
# Author: Paige Smallman, 2026
################################################################
run_permanova_richness <- function(ps_object, nperm = 2999, seed = 123, results_file = "stats_results.txt") {
  suppressPackageStartupMessages({
    require(phyloseq)
    require(vegan)
    require(dplyr)
    require(tibble)
  })

  if (!is.null(seed)) set.seed(seed)
  
  # 1. Setup Metadata
  meta <- as(sample_data(ps_object), "data.frame")
  if (!"Sample" %in% colnames(meta)) {
    meta <- meta %>% tibble::rownames_to_column("Sample")
  }
  
  # Ensure factors are correctly typed
  if ("Year" %in% colnames(meta)) meta$year <- meta$Year
  meta$year   <- as.factor(meta$year)
  meta$season <- as.factor(meta$season)
  meta$local  <- as.factor(meta$local)
  meta$site   <- as.factor(meta$site)
  
  # 2. Calculate Alpha Diversity (Observed Richness)
  # We extract richness first to ensure our distance matrix is univariate
  richness_df <- phyloseq::estimate_richness(ps_object, measures = "Observed")
  richness_df$Sample <- rownames(richness_df)
  
  # Merge richness with metadata for testing
  analysis_df <- left_join(meta, richness_df, by = "Sample")
  
  # 3. Create Euclidean Distance Matrix of RICHNESS
  # This is the "dist" used for PERMANOVA on a single variable
  dist_rich <- dist(analysis_df$Observed, method = "euclidean")
  
  output_lines <- c("### INVERTEBRATE RICHNESS STATISTICAL ANALYSIS ###", "")
  
  # 4. Add ENSO column if date_col exists
  if ("date_col" %in% colnames(analysis_df)) {
    analysis_df <- analysis_df %>%
      mutate(
        Date = as.Date(date_col, format = "%d-%m-%Y"),
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
    analysis_df$ENSO <- as.factor(analysis_df$ENSO)
  }
  
  # 5. Shapiro-Wilk Normality Tests
  output_lines <- c(output_lines, "--- Shapiro-Wilk: Normality of Richness by Year/Location ---")
  for (loc in unique(analysis_df$local)) {
    for (yr in unique(analysis_df$year)) {
      df_sub <- analysis_df[analysis_df$local == loc & analysis_df$year == yr, ]
      if (nrow(df_sub) >= 3) {
        sw <- shapiro.test(df_sub$Observed)
        output_lines <- c(output_lines, paste0(loc, " (", yr, "): W=", round(sw$statistic, 4), " p=", round(sw$p.value, 5)))
      }
    }
  }
  
  # 6. Global PERMANOVA (Interactions + Strata)
  # local * year * season = main effects AND all possible interactions
  # strata = site accounts for your multiple sites/replications
  output_lines <- c(output_lines, "", "--- Global PERMANOVA (Richness ~ Local * Year * Season) ---")
  
  # Use tryCatch to prevent script failure if permutations are impossible
  adonis_global <- tryCatch({
    vegan::adonis2(dist_rich ~ local * year * season, 
                   data = analysis_df, 
                   permutations = nperm, 
                   strata = analysis_df$site,
                   by = "terms")
  }, error = function(e) paste("Error in PERMANOVA:", e$message))
  
  output_lines <- c(output_lines, capture.output(print(adonis_global)))
  
  # 7. PERMDISP (Homogeneity of Multivariate Dispersion)
  # Essential because 2021 failed normality; we must check if variance is equal
  output_lines <- c(output_lines, "", "--- PERMDISP (Testing Homogeneity of Variance by Location) ---")
  disp_loc <- vegan::betadisper(dist_rich, analysis_df$local)
  perm_disp <- vegan::permutest(disp_loc, permutations = nperm)
  output_lines <- c(output_lines, capture.output(print(perm_disp)))
  
  # 7a. PERMANOVA for ENSO
  if ("ENSO" %in% colnames(analysis_df)) {
    output_lines <- c(output_lines, "", "--- PERMANOVA (Richness ~ ENSO * Local) ---")
    adonis_enso <- tryCatch({
      vegan::adonis2(dist_rich ~ ENSO * local, 
                     data = analysis_df, 
                     permutations = nperm, 
                     strata = analysis_df$site,
                     by = "terms")
    }, error = function(e) paste("Error in ENSO PERMANOVA:", e$message))
    output_lines <- c(output_lines, capture.output(print(adonis_enso)))
    
    # Pairwise post-hoc for ENSO within each location
    output_lines <- c(output_lines, "", "--- PAIRWISE POST-HOC: ENSO within Locations ---")
    if (!requireNamespace("pairwiseAdonis", quietly = TRUE)) {
      stop("The pairwiseAdonis package is required for pairwise post-hoc PERMANOVA.")
    }
    for (loc in unique(analysis_df$local)) {
      output_lines <- c(output_lines, paste0("\n>> Testing ENSO differences within Location: ", loc))
      sub_df <- analysis_df[analysis_df$local == loc, ]
      if (length(unique(sub_df$ENSO)) > 1 && nrow(sub_df) > 2) {
        sub_dist <- dist(sub_df$Observed, method = "euclidean")
        pair_enso <- pairwiseAdonis::pairwise.adonis2(sub_dist ~ ENSO, data = sub_df, strata = sub_df$site, nperm = nperm)
        pair_enso$parent_comparison <- NULL
        output_lines <- c(output_lines, capture.output(print(pair_enso)))
        # Note: These pairwise results are printed to the file but not returned in the list object yet.
      }
    }
  }
  
  # 7b. Pairwise Post-hoc PERMANOVA if global is significant for year or local:year
  pairwise_results <- list()
  sig_year <- FALSE
  sig_local_year <- FALSE
  if (inherits(adonis_global, "anova")) {
    # Check significance for 'year' and 'local:year'
    p_year <- NA
    p_local_year <- NA
    if ("year" %in% rownames(adonis_global)) p_year <- adonis_global["year", "Pr(>F)"]
    if ("local:year" %in% rownames(adonis_global)) p_local_year <- adonis_global["local:year", "Pr(>F)"]
    sig_year <- !is.na(p_year) && p_year < 0.05
    sig_local_year <- !is.na(p_local_year) && p_local_year < 0.05
  }
  
  if (sig_year || sig_local_year) {
    output_lines <- c(output_lines, "", "--- PAIRWISE POST-HOC: Years within Locations ---")
    if (!requireNamespace("pairwiseAdonis", quietly = TRUE)) {
      stop("The pairwiseAdonis package is required for pairwise post-hoc PERMANOVA.")
    }
    for (loc in unique(analysis_df$local)) {
      output_lines <- c(output_lines, paste0("\n>> Testing Year differences within Location: ", loc))
      sub_df <- analysis_df[analysis_df$local == loc, ]
      if (length(unique(sub_df$year)) > 1 && nrow(sub_df) > 2) {
        sub_dist <- dist(sub_df$Observed, method = "euclidean")
        pair_year <- pairwiseAdonis::pairwise.adonis2(sub_dist ~ year, data = sub_df, strata = sub_df$site, nperm = nperm)
        pair_year$parent_comparison <- NULL
        output_lines <- c(output_lines, capture.output(print(pair_year)))
        pairwise_results[[as.character(loc)]] <- pair_year
      }
    }
  }
  
  print(adonis_global)
  
  # 8. Final Save and Return (includes pairwise results)
  if (!is.null(results_file) && results_file != "") {
    writeLines(output_lines, results_file)
    message("Results saved to: ", results_file)
  }
  
  # Return results as a list for plotting or further post-hoc tests
  return(list(
    data = analysis_df,
    permanova = adonis_global,
    dispersion = perm_disp,
    pairwise = pairwise_results
  ))
}
################################################################
# (COI-H) RAREFACTION PLOTS - DNA metabarcoding visualization
################################################################
# Purpose: Generate per-sample rarefaction curves from a phyloseq object.
#
# Defaults in this function are set to:
#   - remove control samples using sample_data `sample_type`
#   - use raw ASV counts
#   - save a single PNG output
#
# Author: Paige Smallman, 2026
################################################################

plot_rarefaction_curves <- function(ps_object,
																		figures_dir,
																		project_id,
																		control_col = "sample_type",
																		neg_controls = c("extraction control", "field control", "pcr control"),
												min_step = 50,
												accum_permutations = 1000,
												accum_group_col = "sample_type") {
	suppressPackageStartupMessages({
		require(phyloseq)
		require(ggplot2)
		require(vegan)
		require(dplyr)
	})

	if (!is(ps_object, "phyloseq")) {
		stop("`ps_object` must be a phyloseq object.")
	}

	if (!dir.exists(figures_dir)) {
		dir.create(figures_dir, recursive = TRUE)
	}

	# Fixed palette requested for season plots
	season_palette <- c(
		"Upwelling" = "#56B4E9",
		"Non-Upwelling" = "#E69F00"
	)

	# --- 1. Remove control samples ---
	ps_use <- ps_object
	if (control_col %in% colnames(sample_data(ps_use))) {
		ctrl_values <- tolower(trimws(as.character(sample_data(ps_use)[[control_col]])))
		neg_controls_std <- tolower(trimws(neg_controls))
		keep_samples <- !(ctrl_values %in% neg_controls_std)
		ps_use <- prune_samples(keep_samples, ps_use)
	} else {
		warning("Control column '", control_col, "' not found in sample_data; no controls were removed.")
	}

	# Keep only informative samples/taxa
	ps_use <- prune_samples(sample_sums(ps_use) > 0, ps_use)
	ps_use <- prune_taxa(taxa_sums(ps_use) > 0, ps_use)

	if (nsamples(ps_use) == 0) {
		stop("No samples remain after filtering controls and zero-read samples.")
	}

	# --- 2. Build count matrix (raw ASV counts) ---
	otu_mat <- as(otu_table(ps_use), "matrix")
	if (taxa_are_rows(ps_use)) {
		otu_mat <- t(otu_mat)
	}
	storage.mode(otu_mat) <- "numeric"
	sample_ids <- sample_names(ps_use)

	lib_sizes <- rowSums(otu_mat)
	min_depth <- min(lib_sizes)
	if (min_depth < 2) {
		stop("Minimum sample depth is < 2 reads; cannot compute rarefaction curves reliably.")
	}

	# Ensure sample IDs exist for downstream plotting
	if (is.null(rownames(otu_mat)) || any(rownames(otu_mat) == "")) {
		if (length(sample_ids) == nrow(otu_mat)) {
			rownames(otu_mat) <- sample_ids
		} else {
			rownames(otu_mat) <- paste0("Sample_", seq_len(nrow(otu_mat)))
		}
	}

	step_size <- max(1, min(min_step, floor(min_depth / 10)))

	# --- 3. Compute rarefaction curves per sample ---
	rare_list <- vegan::rarecurve(otu_mat, step = step_size, label = FALSE)
	rare_ids <- names(rare_list)
	if (is.null(rare_ids) || length(rare_ids) != length(rare_list) || any(is.na(rare_ids)) || any(rare_ids == "")) {
		rare_ids <- rownames(otu_mat)
	}

	# Map sample -> Year for curve coloring
	meta_df <- as.data.frame(sample_data(ps_use), stringsAsFactors = FALSE)
	if ("Year" %in% colnames(meta_df)) {
		year_map <- as.character(meta_df$Year)
		year_map[is.na(year_map) | trimws(year_map) == ""] <- "Unknown"
	} else {
		year_map <- rep("Unknown", nrow(meta_df))
	}
	if (!is.null(rownames(meta_df)) && !any(rownames(meta_df) == "")) {
		names(year_map) <- rownames(meta_df)
	} else if (length(sample_ids) == length(year_map)) {
		names(year_map) <- sample_ids
	}

	get_group_map <- function(col_name, fallback = "Unknown", normalize_labels = TRUE) {
		if (col_name %in% colnames(meta_df)) {
			vals <- as.character(meta_df[[col_name]])
		} else {
			vals <- rep(fallback, nrow(meta_df))
		}

		# Normalize spacing/case to avoid duplicate-looking groups
		vals[is.na(vals)] <- ""
		vals <- gsub("\\u00A0", " ", vals) # non-breaking spaces
		vals <- gsub("_", " ", vals)
		vals <- gsub("[[:space:]]+", " ", vals)
		vals <- trimws(vals)
		if (normalize_labels) {
			vals <- tolower(vals)
			vals <- tools::toTitleCase(vals)
		}

		# Canonicalize known grouping labels
		if (tolower(col_name) == "season") {
			vals_l <- tolower(vals)
			vals[grepl("^upwelling$", vals_l)] <- "Upwelling"
			vals[grepl("^non[- ]?upwelling$", vals_l)] <- "Non-Upwelling"
		}
		if (tolower(col_name) == "sample_type") {
			vals_l <- tolower(vals)
			vals[grepl("^edna[ -]?sediment$", vals_l)] <- "eDNA Sediment"
			vals[grepl("^edna[ -]?water$", vals_l)] <- "eDNA Water"
		}

		vals[vals == ""] <- fallback

		if (!is.null(rownames(meta_df)) && !any(rownames(meta_df) == "")) {
			names(vals) <- rownames(meta_df)
		} else if (length(sample_ids) == length(vals)) {
			names(vals) <- sample_ids
		}
		vals
	}

	local_map <- get_group_map("local", fallback = "Unknown", normalize_labels = TRUE)
	season_map <- get_group_map("season", fallback = "Unknown", normalize_labels = TRUE)

	rare_df <- bind_rows(lapply(seq_along(rare_list), function(i) {
		y <- rare_list[[i]]
		reads <- suppressWarnings(as.numeric(names(y)))
		if (all(is.na(reads))) {
			reads <- seq_along(y)
		}
		year_i <- unname(year_map[rare_ids[i]])
		if (length(year_i) == 0 || is.na(year_i) || year_i == "") {
			year_i <- "Unknown"
		}
		local_i <- unname(local_map[rare_ids[i]])
		if (length(local_i) == 0 || is.na(local_i) || local_i == "") {
			local_i <- "Unknown"
		}
		season_i <- unname(season_map[rare_ids[i]])
		if (length(season_i) == 0 || is.na(season_i) || season_i == "") {
			season_i <- "Unknown"
		}
		data.frame(
			Sample = rep(rare_ids[i], length(y)),
			Year = rep(year_i, length(y)),
			Local = rep(local_i, length(y)),
			Season = rep(season_i, length(y)),
			Reads = reads,
			Observed_ASVs = as.numeric(y),
			stringsAsFactors = FALSE
		)
	}))
	rare_df$Year <- as.factor(rare_df$Year)
	rare_df$Local <- as.factor(rare_df$Local)
	rare_df$Season <- as.factor(rare_df$Season)

	# --- 4. Plot and save PNG ---
	p_rare <- ggplot(rare_df, aes(x = Reads, y = Observed_ASVs, group = Sample, color = Year)) +
		geom_line(alpha = 0.20, linewidth = 0.7) +
		theme_minimal() +
		labs(
			title = paste(project_id, "- Rarefaction Curves by Year (Per Sample, Controls Removed)"),
			x = "Sequencing Depth (Reads)",
			y = "Observed ASV Richness",
			color = "Year"
		)

	out_png <- file.path(figures_dir, paste0(project_id, "_rarefaction_per_sample.png"))
	ggsave(out_png, p_rare, width = 10, height = 7, dpi = 300)

	# --- 4b. Rarefaction curves colored by Location ---
	p_rare_local <- ggplot(rare_df, aes(x = Reads, y = Observed_ASVs, group = Sample, color = Local)) +
		geom_line(alpha = 0.20, linewidth = 0.7) +
		theme_minimal() +
		labs(
			title = paste(project_id, "- Rarefaction Curves by Location (Per Sample, Controls Removed)"),
			x = "Sequencing Depth (Reads)",
			y = "Observed ASV Richness",
			color = "Location"
		)
	out_local_png <- file.path(figures_dir, paste0(project_id, "_rarefaction_by_local.png"))
	ggsave(out_local_png, p_rare_local, width = 10, height = 7, dpi = 300)

	# --- 4c. Rarefaction curves colored by Season ---
	p_rare_season <- ggplot(rare_df, aes(x = Reads, y = Observed_ASVs, group = Sample, color = Season)) +
		geom_line(alpha = 0.20, linewidth = 0.7) +
		scale_color_manual(values = season_palette, drop = FALSE) +
		theme_minimal() +
		labs(
			title = paste(project_id, "- Rarefaction Curves by Season (Per Sample, Controls Removed)"),
			x = "Sequencing Depth (Reads)",
			y = "Observed ASV Richness",
			color = "Season"
		)
	out_season_png <- file.path(figures_dir, paste0(project_id, "_rarefaction_by_season.png"))
	ggsave(out_season_png, p_rare_season, width = 10, height = 7, dpi = 300)

	# --- 5. Sample accumulation curve (samples on x, ASV richness on y) ---
	otu_pa <- otu_mat
	otu_pa[otu_pa > 0] <- 1

	acc <- vegan::specaccum(otu_pa, method = "random", permutations = accum_permutations)
	acc_df <- data.frame(
		Samples = acc$sites,
		Richness = acc$richness,
		SD = acc$sd,
		stringsAsFactors = FALSE
	)
	acc_df$Lower <- pmax(0, acc_df$Richness - 1.96 * acc_df$SD)
	acc_df$Upper <- acc_df$Richness + 1.96 * acc_df$SD

	p_acc <- ggplot(acc_df, aes(x = Samples, y = Richness)) +
		geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#A6CEE3", alpha = 0.35) +
		geom_line(color = "#1F78B4", linewidth = 0.9) +
		theme_minimal() +
		labs(
			title = paste(project_id, "- Sample Accumulation Curve (Controls Removed)"),
			x = "Number of Samples",
			y = "Accumulated ASV Richness"
		)

	out_acc_png <- file.path(figures_dir, paste0(project_id, "_sample_accumulation_curve.png"))
	ggsave(out_acc_png, p_acc, width = 10, height = 7, dpi = 300)

	# --- 6. Grouped sample accumulation curves ---
	plot_grouped_accumulation <- function(group_col) {
		normalize_labels <- !(tolower(group_col) %in% c("year"))
		acc_group_map <- get_group_map(group_col, fallback = "Unknown", normalize_labels = normalize_labels)
		acc_group_vec <- unname(acc_group_map[rownames(otu_pa)])
		acc_group_vec[is.na(acc_group_vec) | trimws(acc_group_vec) == ""] <- "Unknown"

		group_levels <- sort(unique(acc_group_vec))
		acc_group_df <- bind_rows(lapply(group_levels, function(g) {
			idx <- which(acc_group_vec == g)
			if (length(idx) < 2) {
				return(NULL)
			}
			otu_g <- otu_pa[idx, , drop = FALSE]
			acc_g <- vegan::specaccum(otu_g, method = "random", permutations = accum_permutations)
			data.frame(
				Group = g,
				Samples = acc_g$sites,
				Richness = acc_g$richness,
				SD = acc_g$sd,
				stringsAsFactors = FALSE
			)
		}))

		if (nrow(acc_group_df) > 0) {
			acc_group_df$Lower <- pmax(0, acc_group_df$Richness - 1.96 * acc_group_df$SD)
			acc_group_df$Upper <- acc_group_df$Richness + 1.96 * acc_group_df$SD

			p_acc_group <- ggplot(acc_group_df, aes(x = Samples, y = Richness, color = Group, fill = Group, group = Group)) +
				geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, linewidth = 0, show.legend = FALSE) +
				geom_line(linewidth = 1.0) +
				theme_minimal() +
				labs(
					title = paste(project_id, "- Sample Accumulation by", group_col, "(Controls Removed)"),
					x = "Number of Samples",
					y = "Accumulated ASV Richness",
					color = group_col
				)

			if (tolower(group_col) == "season") {
				p_acc_group <- p_acc_group +
					scale_color_manual(values = season_palette, drop = FALSE) +
					scale_fill_manual(values = season_palette, drop = FALSE)
			}

			safe_group_col <- gsub("[^A-Za-z0-9]+", "_", group_col)
			out_acc_group_png <- file.path(figures_dir, paste0(project_id, "_sample_accumulation_by_", safe_group_col, ".png"))
			ggsave(out_acc_group_png, p_acc_group, width = 10, height = 7, dpi = 300)
			message("Sample accumulation by ", group_col, " saved: ", out_acc_group_png)
		} else {
			warning("No grouped sample accumulation curves were generated for '", group_col, "' (need at least 2 samples per group).")
		}
	}

	group_cols_to_plot <- unique(c(accum_group_col, "local", "season", "Year"))
	for (gcol in group_cols_to_plot) {
		plot_grouped_accumulation(gcol)
	}

	message("Rarefaction plot saved: ", out_png)
	message("Rarefaction plot by location saved: ", out_local_png)
	message("Rarefaction plot by season saved: ", out_season_png)
	message("Sample accumulation plot saved: ", out_acc_png)
	invisible(p_rare)
}
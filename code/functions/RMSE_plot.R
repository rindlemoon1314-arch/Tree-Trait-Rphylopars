plot_family_trait_rmse_phylogeny <- function(
    tree = here::here("Tree Trait data/no_monocots_tree.nwk"),
    trait_data = here::here("Tree Trait data/TRY_trait_data_cleaned.csv"),
    taxonomy = here::here("V.PhyloMaker2/data/tips.info.TPL.rda"),
    phy_results_dir = c(
      here::here("Results/kfold_cv/single_traits_PHY"),
      here::here("Results/kfold_cv/single_traits_all8_PHY")
    ),
    trait_names = c(
      "Wood_density", "Specific_leaf_area", "Seed_dry_mass",
      "Leaf_P_per_mass", "Stem_conduit_diameter", "Tree_height",
      "Root_depth", "Bark_thickness", "Leaf_N_per_mass",
      "Leaf_thickness", "Leaf_Vcmax_per_dry_mass",
      "Stomatal_conductance", "Leaf_area", "Leaf_K_per_mass",
      "Leaf_density", "Stem_diameter", "Crown_height",
      "Crown_diameter"
    ),
    species_selection_traits = trait_names,
    model = "lambda",
    configuration = "TTF",
    folds = c(1L, 3L, 5L, 6L, 8L, 10L),
    truth_column = "transformed_truth",
    output_file = NULL,
    rmse_file = NULL,
    tip_metadata_file = NULL,
    family_palette_file = NULL,
    plot_title = "Family-level cross-validation error across selected traits",
    width = 8.25,
    height = 10.5,
    dpi = 300,
    base_pointsize = 10.5,
    family_edge_width = 0.42,
    mixed_edge_width = 0.18,
    mixed_edge_colour = "#BEBEBE",
    family_separator_colour = "#4D4D4D",
    family_separator_width = 0.28,
    ring_gap_fraction = 0.0022,
    ring_width_fraction = 0.0165,
    tree_to_ring_gap_fraction = 0.007,
    radial_compression_power = 1.70,
    family_ring_width_fraction = 0.012,
    largest_family_labels = 20L,
    family_key_groups = 8L,
    show_complete_family_key = FALSE,
    show_tip_labels = FALSE,
    tip_label_cex = 0.08) {

  # Use a fixed Word-compatible figure width and prevent the figure from
  # exceeding the available page height.
  width <- 8.25
  height <- min(as.numeric(height), 11.25)
  if (!is.finite(height) || height <= 0) {
    stop("height must be a positive finite number.", call. = FALSE)
  }

  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required.", call. = FALSE)
  }

  normalise_species <- function(x) {
    x <- trimws(as.character(x))
    gsub("\\s+", "_", x)
  }

  normalise_trait <- function(x) {
    x <- trimws(as.character(x))
    x <- gsub("[._/]+", " ", x)
    x <- gsub("\\s+", " ", x)
    tolower(x)
  }

  pretty_trait <- function(x) {
    x <- gsub("_", " ", x, fixed = TRUE)
    x <- gsub("lambda", "lambda", x, fixed = TRUE)
    x
  }

  read_tree_object <- function(x) {
    if (inherits(x, "phylo")) return(x)
    if (!is.character(x) || length(x) != 1L || !file.exists(x)) {
      stop("tree must be a phylo object or an existing file.", call. = FALSE)
    }
    ext <- tolower(tools::file_ext(x))
    if (ext == "rds") return(readRDS(x))
    if (ext %in% c("rda", "rdata")) {
      e <- new.env(parent = emptyenv())
      loaded <- load(x, envir = e)
      candidates <- loaded[vapply(loaded, function(n) inherits(e[[n]], "phylo"), logical(1))]
      if (!length(candidates)) stop("No phylo object was found in tree file.", call. = FALSE)
      return(e[[candidates[[1L]]]])
    }
    ape::read.tree(x)
  }

  read_trait_data <- function(x) {
    if (is.data.frame(x)) return(x)
    if (!is.character(x) || length(x) != 1L || !file.exists(x)) {
      stop("trait_data must be a data frame or an existing CSV/RDS file.", call. = FALSE)
    }
    if (tolower(tools::file_ext(x)) == "rds") return(readRDS(x))
    utils::read.csv(x, stringsAsFactors = FALSE, check.names = FALSE)
  }

  read_taxonomy <- function(x) {
    if (is.data.frame(x)) return(x)
    if (!is.character(x) || length(x) != 1L || !file.exists(x)) {
      stop("taxonomy must be a data frame or an existing file.", call. = FALSE)
    }
    ext <- tolower(tools::file_ext(x))
    if (ext == "csv") return(utils::read.csv(x, stringsAsFactors = FALSE, check.names = FALSE))
    if (ext == "rds") return(readRDS(x))
    e <- new.env(parent = emptyenv())
    loaded <- load(x, envir = e)
    candidates <- loaded[vapply(loaded, function(n) {
      z <- e[[n]]
      is.data.frame(z) && all(c("species", "family") %in% names(z))
    }, logical(1))]
    if (!length(candidates)) stop("No taxonomy table with species and family was found.", call. = FALSE)
    e[[candidates[[1L]]]]
  }

  make_family_palette <- function(families) {
    families <- sort(unique(families))
    hues <- (seq_along(families) * 137.507764) %% 360
    stats::setNames(grDevices::hcl(h = hues, c = 70, l = 52, fixup = TRUE), families)
  }

  assign_node_families <- function(current_tree, tip_families) {
    ntip <- ape::Ntip(current_tree)
    total <- ntip + current_tree$Nnode
    node_family <- rep(NA_character_, total)
    node_family[seq_len(ntip)] <- tip_families
    unresolved <- seq.int(ntip + 1L, total)
    while (length(unresolved)) {
      resolved <- integer()
      for (node in unresolved) {
        children <- current_tree$edge[current_tree$edge[, 1L] == node, 2L]
        child_family <- node_family[children]
        if (length(children) && all(!is.na(child_family))) {
          node_family[node] <- if (length(unique(child_family)) == 1L) child_family[1L] else "Mixed"
          resolved <- c(resolved, node)
        }
      }
      if (!length(resolved)) break
      unresolved <- setdiff(unresolved, resolved)
    }
    node_family[is.na(node_family)] <- "Mixed"
    node_family
  }

  rmse_class <- function(x) {
    result <- rep("No data", length(x))
    result[!is.na(x) & x >= 0 & x <= 0.25] <- "0-0.25"
    result[!is.na(x) & x > 0.25 & x <= 0.50] <- "0.25-0.50"
    result[!is.na(x) & x > 0.50 & x <= 0.75] <- "0.50-0.75"
    result[!is.na(x) & x > 0.75 & x <= 1.00] <- "0.75-1.00"
    result[!is.na(x) & x > 1.00] <- ">1.00"
    factor(result, levels = c("No data", "0-0.25", "0.25-0.50", "0.50-0.75", "0.75-1.00", ">1.00"))
  }

  below_one_palette <- grDevices::colorRampPalette(
    c("#EFF6FF", "#BFD7EA", "#6BAED6", "#08519C")
  )(256L)
  above_one_colour <- "#99000D"

  map_rmse_colour <- function(x) {
    answer <- rep("#FFFFFF", length(x))
    below <- is.finite(x) & x <= 1
    above <- is.finite(x) & x > 1
    if (any(below)) {
      index <- 1L + floor(pmin(pmax(x[below], 0), 1) * 255)
      answer[below] <- below_one_palette[index]
    }
    answer[above] <- above_one_colour
    answer
  }

  full_tree <- read_tree_object(tree)
  full_tree$tip.label <- normalise_species(full_tree$tip.label)
  raw_traits <- read_trait_data(trait_data)
  required_trait_columns <- c("accepted_bin", "trait", "value")
  if (!all(required_trait_columns %in% names(raw_traits))) {
    stop("Long trait_data must contain accepted_bin, trait, and value.", call. = FALSE)
  }

  trait_names <- unique(as.character(trait_names))
  species_selection_traits <- unique(as.character(species_selection_traits))
  selected_keys <- normalise_trait(species_selection_traits)
  observed_rows <- normalise_trait(raw_traits$trait) %in% selected_keys & !is.na(raw_traits$value)
  observed_species <- unique(normalise_species(raw_traits$accepted_bin[observed_rows]))
  observed_species <- observed_species[!is.na(observed_species) & nzchar(observed_species)]
  matched_species <- intersect(full_tree$tip.label, observed_species)
  if (length(matched_species) < 2L) stop("Fewer than two observed species matched the tree.", call. = FALSE)

  selected_tree <- ape::keep.tip(full_tree, matched_species)
  selected_tree <- ape::reorder.phylo(selected_tree, "postorder")

  tax <- read_taxonomy(taxonomy)
  tax$species <- normalise_species(tax$species)
  tax$family <- trimws(as.character(tax$family))
  tax <- tax[!is.na(tax$species) & nzchar(tax$species) & !is.na(tax$family) & nzchar(tax$family), c("species", "family")]
  tax <- tax[!duplicated(tax$species), ]
  family_by_species <- stats::setNames(tax$family, tax$species)
  tax_genus <- sub("_.*$", "", tax$species)
  genus_table <- unique(data.frame(
    genus = tax_genus,
    family = tax$family,
    stringsAsFactors = FALSE
  ))
  genus_table <- genus_table[!duplicated(genus_table$genus), , drop = FALSE]
  family_by_genus <- stats::setNames(genus_table$family, genus_table$genus)
  tip_family <- unname(family_by_species[selected_tree$tip.label])
  missing_tip_family <- is.na(tip_family) | !nzchar(tip_family)
  tip_genus <- sub("_.*$", "", selected_tree$tip.label)
  tip_family[missing_tip_family] <- unname(
    family_by_genus[tip_genus[missing_tip_family]]
  )
  tip_family[is.na(tip_family) | !nzchar(tip_family)] <- "Unknown"

  fold_names <- paste0("fold_", as.integer(folds))
  prediction_parts <- list()
  part_index <- 0L

  for (trait in trait_names) {
    result_basename <- paste0(
      "PHY_kfoldcv_", trait, "_", model, "_", configuration, ".rds"
    )
    result_candidates <- file.path(phy_results_dir, result_basename)
    result_candidates <- result_candidates[file.exists(result_candidates)]
    if (!length(result_candidates)) {
      warning("Missing PHY file: ", result_basename, call. = FALSE)
      next
    }
    result_file <- result_candidates[[1L]]
    result <- readRDS(result_file)
    status <- attr(result, "status")
    assessments <- attr(result, "assessment")
    if (!is.data.frame(status) || !is.list(assessments)) {
      warning("Invalid PHY metadata: ", basename(result_file), call. = FALSE)
      next
    }
    for (fold_name in fold_names) {
      fold_index <- match(fold_name, as.character(status$fold))
      if (is.na(fold_index) || status$status[fold_index] != "success" || is.null(result[[fold_index]])) next
      fold_result <- result[[fold_index]]
      assessment <- assessments[[fold_index]]
      if (!is.data.frame(assessment) || !all(c("species", "trait", truth_column) %in% names(assessment))) next
      available_trait <- intersect(colnames(fold_result$anc_recon), colnames(fold_result$anc_var))
      if (!(trait %in% available_trait)) next
      assessment <- assessment[
        assessment$trait == trait & is.finite(assessment[[truth_column]]),
        c("species", "trait", truth_column), drop = FALSE
      ]
      if (!nrow(assessment)) next
      names(assessment)[3L] <- "truth"
      assessment$species <- normalise_species(assessment$species)
      truth <- stats::aggregate(truth ~ species + trait, assessment, mean)
      rows <- match(truth$species, rownames(fold_result$anc_recon))
      cols <- match(truth$trait, colnames(fold_result$anc_recon))
      valid <- !is.na(rows) & !is.na(cols)
      if (!any(valid)) next
      truth <- truth[valid, , drop = FALSE]
      prediction <- fold_result$anc_recon[cbind(rows[valid], cols[valid])]
      part_index <- part_index + 1L
      prediction_parts[[part_index]] <- data.frame(
        fold = fold_name,
        species = truth$species,
        trait = trait,
        truth = truth$truth,
        prediction = as.numeric(prediction),
        squared_error = (as.numeric(prediction) - truth$truth)^2,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(prediction_parts)) stop("No valid CV predictions were extracted.", call. = FALSE)
  predictions <- do.call(rbind, prediction_parts)
  predictions$family <- unname(family_by_species[predictions$species])
  missing_prediction_family <- is.na(predictions$family) |
    !nzchar(predictions$family)
  prediction_genus <- sub("_.*$", "", predictions$species)
  predictions$family[missing_prediction_family] <- unname(
    family_by_genus[prediction_genus[missing_prediction_family]]
  )
  predictions <- predictions[!is.na(predictions$family) & nzchar(predictions$family), , drop = FALSE]

  rmse_table <- stats::aggregate(
    squared_error ~ family + trait,
    predictions,
    function(x) sqrt(mean(x))
  )
  names(rmse_table)[names(rmse_table) == "squared_error"] <- "rmse"
  count_table <- stats::aggregate(
    squared_error ~ family + trait,
    predictions,
    length
  )
  names(count_table)[names(count_table) == "squared_error"] <- "n_predictions"
  rmse_table <- merge(
    rmse_table,
    count_table,
    by = c("family", "trait"),
    all = TRUE,
    sort = FALSE
  )
  rmse_table$rmse_class <- rmse_class(rmse_table$rmse)

  all_families <- sort(unique(tip_family))
  grid <- expand.grid(family = all_families, trait = trait_names, stringsAsFactors = FALSE)
  grid_key <- paste(grid$family, grid$trait, sep = "\r")
  rmse_key <- paste(rmse_table$family, rmse_table$trait, sep = "\r")
  matched <- match(grid_key, rmse_key)
  grid$rmse <- rmse_table$rmse[matched]
  grid$n_predictions <- rmse_table$n_predictions[matched]
  grid$rmse_class <- rmse_class(grid$rmse)
  rmse_table_complete <- grid

  rmse_table_complete$rmse_colour <- map_rmse_colour(rmse_table_complete$rmse)

  family_palette <- make_family_palette(setdiff(all_families, "Unknown"))
  node_family <- assign_node_families(selected_tree, tip_family)
  edge_family <- node_family[selected_tree$edge[, 2L]]
  edge_colour <- rep(mixed_edge_colour, nrow(selected_tree$edge))
  family_edge <- edge_family %in% names(family_palette)
  edge_colour[family_edge] <- unname(family_palette[edge_family[family_edge]])
  edge_width <- ifelse(family_edge, family_edge_width, mixed_edge_width)

  open_device <- function(path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    ext <- tolower(tools::file_ext(path))
    if (ext == "pdf") {
      grDevices::cairo_pdf(path, width = width, height = height,
                           pointsize = base_pointsize)
    } else if (ext == "png") {
      grDevices::png(path, width = width, height = height, units = "in",
                     res = dpi, type = "cairo", pointsize = base_pointsize)
    } else stop("output_file must end in .pdf or .png.", call. = FALSE)
  }

  annular_sector <- function(a0, a1, r0, r1, colour, border = NA) {
    n <- max(3L, ceiling(abs(a1 - a0) * 60))
    theta <- seq(a0, a1, length.out = n)
    graphics::polygon(
      c(r0 * cos(theta), rev(r1 * cos(theta))),
      c(r0 * sin(theta), rev(r1 * sin(theta))),
      col = colour, border = border
    )
  }

  draw_plot <- function() {
    layout(
      matrix(c(1, 2), ncol = 1),
      heights = if (isTRUE(show_complete_family_key)) c(6.65, 6.60) else c(6.65, 0.95)
    )
    on.exit(layout(1), add = TRUE)

    if (!is.numeric(radial_compression_power) ||
        length(radial_compression_power) != 1L ||
        !is.finite(radial_compression_power) ||
        radial_compression_power < 1) {
      stop("radial_compression_power must be a finite number >= 1.", call. = FALSE)
    }
    plot_tree <- selected_tree
    if (radial_compression_power > 1) {
      original_depth <- ape::node.depth.edgelength(plot_tree)
      maximum_depth <- max(original_depth, na.rm = TRUE)
      compressed_depth <- maximum_depth *
        (pmax(original_depth, 0) / maximum_depth)^radial_compression_power
      parent_depth <- compressed_depth[plot_tree$edge[, 1L]]
      child_depth <- compressed_depth[plot_tree$edge[, 2L]]
      plot_tree$edge.length <- pmax(child_depth - parent_depth, .Machine$double.eps)
    }
    depth <- ape::node.depth.edgelength(plot_tree)
    tree_radius <- max(depth[seq_len(ape::Ntip(selected_tree))], na.rm = TRUE)
    ring_gap <- tree_radius * ring_gap_fraction
    ring_width <- tree_radius * ring_width_fraction
    first_ring <- tree_radius * (1 + tree_to_ring_gap_fraction)
    trait_outer_radius <- first_ring + length(trait_names) * (ring_width + ring_gap)
    family_ring_gap <- tree_radius * 0.006
    family_ring_width <- tree_radius * family_ring_width_fraction
    family_ring_inner <- trait_outer_radius + family_ring_gap
    family_ring_outer <- family_ring_inner + family_ring_width
    outer_radius <- family_ring_outer + tree_radius * 0.055

    has_plot_title <- !is.null(plot_title) &&
      length(plot_title) == 1L &&
      nzchar(trimws(plot_title))
    graphics::par(
      mar = c(0.05, 0.05, if (has_plot_title) 1.65 else 0.15, 0.05),
      bg = "white"
    )
    ape::plot.phylo(
      plot_tree, type = "fan", show.tip.label = show_tip_labels,
      cex = tip_label_cex, edge.color = edge_colour, edge.width = edge_width,
      no.margin = FALSE, label.offset = 0,
      x.lim = c(-outer_radius * 1.02, outer_radius * 1.02),
      y.lim = c(-outer_radius * 1.02, outer_radius * 1.02)
    )
    plot_environment <- get(".PlotPhyloEnv", envir = asNamespace("ape"))
    last_plot <- get("last_plot.phylo", envir = plot_environment)
    ntip <- ape::Ntip(selected_tree)
    root_node <- ntip + 1L
    cx <- last_plot$xx[root_node]
    cy <- last_plot$yy[root_node]
    tip_angle <- atan2(last_plot$yy[seq_len(ntip)] - cy, last_plot$xx[seq_len(ntip)] - cx)
    tip_angle[tip_angle < 0] <- tip_angle[tip_angle < 0] + 2 * pi
    order_index <- order(tip_angle)
    angles <- tip_angle[order_index]
    ordered_family <- tip_family[order_index]

    run_id <- cumsum(c(TRUE, ordered_family[-1L] != ordered_family[-length(ordered_family)]))
    runs <- split(seq_along(angles), run_id)
    run_bounds <- lapply(runs, function(pos) {
      first <- min(pos); last <- max(pos)
      previous_angle <- if (first == 1L) angles[length(angles)] - 2 * pi else angles[first - 1L]
      next_angle <- if (last == length(angles)) angles[1L] + 2 * pi else angles[last + 1L]
      c(start = (previous_angle + angles[first]) / 2, end = (angles[last] + next_angle) / 2)
    })

    rmse_lookup <- stats::setNames(
      rmse_table_complete$rmse_colour,
      paste(rmse_table_complete$family, rmse_table_complete$trait, sep = "\r")
    )
    for (run_number in seq_along(runs)) {
      positions <- runs[[run_number]]
      family <- ordered_family[positions[1L]]
      bounds <- run_bounds[[run_number]]
      for (trait_index in seq_along(trait_names)) {
        r0 <- first_ring + (trait_index - 1L) * (ring_width + ring_gap)
        r1 <- r0 + ring_width
        key <- paste(family, trait_names[trait_index], sep = "\r")
        colour <- rmse_lookup[[key]]
        if (is.null(colour) || is.na(colour)) colour <- "#FFFFFF"
        annular_sector(bounds["start"], bounds["end"], r0, r1, colour)
      }
      family_colour <- if (family %in% names(family_palette)) {
        family_palette[[family]]
      } else {
        "#D9D9D9"
      }
      annular_sector(
        bounds["start"], bounds["end"], family_ring_inner,
        family_ring_outer, family_colour
      )
    }

    run_family <- vapply(runs, function(pos) ordered_family[pos[1L]], character(1))
    run_mid <- vapply(run_bounds, function(z) mean(z), numeric(1))
    family_counts <- sort(table(tip_family), decreasing = TRUE)
    largest_families <- names(head(family_counts, largest_family_labels))
    labelled_once <- character()
    for (i in seq_along(run_family)) {
      family <- run_family[i]
      if (!family %in% largest_families || family %in% labelled_once) next
      labelled_once <- c(labelled_once, family)
      angle <- run_mid[i]
      label_radius <- family_ring_outer + tree_radius * 0.022
      x0 <- family_ring_outer * cos(angle)
      y0 <- family_ring_outer * sin(angle)
      x1 <- label_radius * cos(angle)
      y1 <- label_radius * sin(angle)
      graphics::segments(x0, y0, x1, y1, col = "#555555", lwd = 0.35)
      degrees <- angle * 180 / pi
      flip <- degrees > 90 && degrees < 270
      graphics::text(
        x1, y1, family,
        srt = if (flip) degrees + 180 else degrees,
        adj = if (flip) c(1.03, 0.5) else c(-0.03, 0.5),
        cex = 0.68, col = "#333333"
      )
    }

    ordered_unique_families <- unique(run_family)
    key_group <- cut(
      seq_along(ordered_unique_families),
      breaks = family_key_groups,
      labels = FALSE,
      include.lowest = TRUE
    )
    family_key <- split(ordered_unique_families, key_group)
    group_for_run <- match(run_family, ordered_unique_families)
    group_for_run <- key_group[group_for_run]
    if (isTRUE(show_complete_family_key)) {
      for (group_index in seq_along(family_key)) {
        run_positions <- which(group_for_run == group_index)
        if (!length(run_positions)) next
        angle <- mean(range(run_mid[run_positions]))
        marker_radius <- family_ring_outer + tree_radius * 0.010
        graphics::segments(
          family_ring_outer * cos(angle), family_ring_outer * sin(angle),
          marker_radius * cos(angle), marker_radius * sin(angle),
          col = "#444444", lwd = 0.38
        )
        graphics::text(
          marker_radius * cos(angle), marker_radius * sin(angle),
          paste0("S", group_index), cex = 0.48, font = 2,
          col = "#111111"
        )
      }
    }

    for (bounds in run_bounds) {
      graphics::segments(
        tree_radius * 0.992 * cos(bounds["start"]),
        tree_radius * 0.992 * sin(bounds["start"]),
        outer_radius * 1.006 * cos(bounds["start"]),
        outer_radius * 1.006 * sin(bounds["start"]),
        col = family_separator_colour, lwd = family_separator_width
      )
    }
    if (has_plot_title) {
      graphics::title(
        main = plot_title,
        cex.main = 1.20,
        font.main = 2,
        line = 0.32
      )
    }

    graphics::par(mar = c(0.12, 0.18, 0.03, 0.18), xpd = NA)
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
    graphics::rect(0.004, 0.025, 0.996, 0.982, col = "white",
                   border = "#B8B8B8", lwd = 0.55)
    legend_header_y <- if (isTRUE(show_complete_family_key)) 0.925 else 0.84
    divider_bottom <- if (isTRUE(show_complete_family_key)) 0.61 else 0.10
    divider_top <- if (isTRUE(show_complete_family_key)) 0.945 else 0.90
    scale_mid_y <- if (isTRUE(show_complete_family_key)) 0.805 else 0.52
    scale_y0 <- if (isTRUE(show_complete_family_key)) 0.780 else 0.43
    scale_y1 <- if (isTRUE(show_complete_family_key)) 0.830 else 0.60
    scale_tick_y <- if (isTRUE(show_complete_family_key)) 0.755 else 0.35

    graphics::text(0.025, legend_header_y, "Outer RMSE rings (inner to outer)",
                   adj = c(0, 0.5), cex = 0.88, font = 2)

    columns <- 2L
    rows <- ceiling(length(trait_names) / columns)
    column_x <- c(0.028, 0.265)
    number_x <- column_x
    label_x <- column_x + 0.035
    trait_y <- if (isTRUE(show_complete_family_key)) {
      seq(0.82, 0.64, length.out = rows)
    } else {
      seq(0.60, 0.20, length.out = rows)
    }
    for (i in seq_along(trait_names)) {
      column <- (i - 1L) %/% rows + 1L
      row <- (i - 1L) %% rows + 1L
      graphics::text(number_x[column], trait_y[row], sprintf("%02d", i),
                     adj = c(0, 0.5), cex = 0.76, col = "#707070")
      graphics::text(label_x[column], trait_y[row], pretty_trait(trait_names[i]),
                     adj = c(0, 0.5), cex = 0.80)
    }

    graphics::segments(0.515, divider_bottom, 0.515, divider_top,
                       col = "#D0D0D0", lwd = 0.55)
    graphics::text(0.545, legend_header_y, "Family-trait RMSE", adj = c(0, 0.5),
                   cex = 0.88, font = 2)
    draw_gradient <- function(x0, x1, y0, y1, palette) {
      breaks <- seq(x0, x1, length.out = length(palette) + 1L)
      for (i in seq_along(palette)) {
        graphics::rect(breaks[i], y0, breaks[i + 1L], y1,
                       col = palette[i], border = NA)
      }
      graphics::rect(x0, y0, x1, y1, border = "#777777", lwd = 0.45)
    }
    graphics::text(0.545, scale_mid_y, "0–1", adj = c(0, 0.5), cex = 0.80)
    draw_gradient(0.595, 0.745, scale_y0, scale_y1, below_one_palette)
    graphics::text(0.595, scale_tick_y, "0", adj = c(0.5, 1), cex = 0.75,
                   col = "#555555")
    graphics::text(0.745, scale_tick_y, "1", adj = c(0.5, 1), cex = 0.75,
                   col = "#555555")
    graphics::text(0.775, scale_mid_y, ">1", adj = c(0, 0.5), cex = 0.80)
    graphics::rect(0.815, scale_y0, 0.925, scale_y1,
                   col = above_one_colour, border = "#777777", lwd = 0.45)
    graphics::text(0.87, scale_tick_y, "All values >1", adj = c(0.5, 1),
                   cex = 0.73, col = "#555555")
    graphics::rect(0.945, scale_y0, 0.970, scale_y1, col = "#FFFFFF",
                   border = "#777777", lwd = 0.45)
    graphics::text(0.9575, scale_tick_y, "No data", adj = c(0.5, 1), cex = 0.70)

    if (isTRUE(show_complete_family_key)) {
      graphics::text(0.025, 0.565,
                     paste0("Family-name sectors (clockwise; S1–S", family_key_groups,
                            " are marked on the outer family ring)"),
                     adj = c(0, 0.5), cex = 0.82, font = 2)
      box_x <- c(0.025, 0.272, 0.519, 0.766)
      box_w <- 0.218

      # Divide the eight family-name sectors with simple rules rather than boxes.
      graphics::segments(
        x0 = c(0.260, 0.507, 0.754),
        y0 = 0.065,
        x1 = c(0.260, 0.507, 0.754),
        y1 = 0.525,
        col = "#C7C7C7", lwd = 0.55
      )
      graphics::segments(
        x0 = 0.025, y0 = 0.297,
        x1 = 0.984, y1 = 0.297,
        col = "#C7C7C7", lwd = 0.55
      )
      for (group_index in seq_along(family_key)) {
        column <- (group_index - 1L) %% 4L + 1L
        row <- (group_index - 1L) %/% 4L + 1L
        x0 <- box_x[column]
        y_top <- if (row == 1L) 0.525 else 0.292
        y_bottom <- if (row == 1L) 0.305 else 0.065
        graphics::text(x0 + 0.010, y_top - 0.022, paste0("S", group_index),
                       adj = c(0, 0.5), cex = 0.80, font = 2)
        names_in_group <- family_key[[group_index]]
        wrapped <- paste(names_in_group, collapse = " · ")
        wrapped <- paste(strwrap(wrapped, width = 42), collapse = "\n")
        graphics::text(x0 + 0.010, y_top - 0.048, wrapped,
                       adj = c(0, 1), cex = 0.47, col = "#252525")
      }
    }
  }

  if (!is.null(output_file)) {
    open_device(output_file)
    draw_plot()
    grDevices::dev.off()
  } else draw_plot()

  tip_metadata <- data.frame(
    species = selected_tree$tip.label,
    family = tip_family,
    stringsAsFactors = FALSE
  )
  palette_table <- data.frame(
    family = names(family_palette),
    colour = unname(family_palette),
    species_count = as.integer(table(factor(tip_family, levels = names(family_palette)))),
    stringsAsFactors = FALSE
  )
  if (!is.null(rmse_file)) {
    dir.create(dirname(rmse_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(rmse_table_complete, rmse_file, row.names = FALSE)
  }
  if (!is.null(tip_metadata_file)) {
    dir.create(dirname(tip_metadata_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(tip_metadata, tip_metadata_file, row.names = FALSE)
  }
  if (!is.null(family_palette_file)) {
    dir.create(dirname(family_palette_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(palette_table, family_palette_file, row.names = FALSE)
  }

  message("Analysis-tree tips: ", format(ape::Ntip(full_tree), big.mark = ","))
  message("Matched and plotted species: ", format(ape::Ntip(selected_tree), big.mark = ","))
  message("Families: ", format(length(unique(tip_family)), big.mark = ","))
  message("CV predictions used: ", format(nrow(predictions), big.mark = ","))
  message("Family-trait combinations with RMSE: ", format(sum(!is.na(rmse_table_complete$rmse)), big.mark = ","))

  invisible(list(
    plot = draw_plot,
    tree = selected_tree,
    predictions = predictions,
    family_trait_rmse = rmse_table_complete,
    tip_metadata = tip_metadata,
    family_palette = palette_table,
    rmse_colours = list(
      below_one = below_one_palette,
      above_one = above_one_colour,
      no_data = "#FFFFFF",
      radial_compression_power = radial_compression_power
    ),
    output_file = output_file
  ))
}

plot_fold_status <- function(
    single_dir = "",
    all8_single_dir = "",
    all8_dir = "",
    high_correlation_dir = "",
    output_file = NULL,
    single_cell_width = 0.78,
    single_cell_height = 0.72,
    multi_cell_width = 0.98,
    multi_cell_height = 0.92,
    figure_width = 8.25,
    figure_height = 11.25,
    trait_text_pt = 9
) {
  required_packages <- c("ggplot2", "patchwork")

  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(missing_packages) > 0L) {
    stop(
      "Required packages are missing: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  result_dirs <- c(
    cluster_single = single_dir,
    all8_single = all8_single_dir,
    all8_multi = all8_dir,
    cluster_multi = high_correlation_dir
  )

  invalid_paths <- (
    is.na(result_dirs) |
      !nzchar(trimws(result_dirs)) |
      !dir.exists(result_dirs)
  )

  if (any(invalid_paths)) {
    stop(
      "These result directories do not exist: ",
      paste(result_dirs[invalid_paths], collapse = ", "),
      call. = FALSE
    )
  }

  validate_number <- function(value, name, lower, upper) {
    if (
      length(value) != 1L ||
      is.na(value) ||
      !is.numeric(value) ||
      value < lower ||
      value > upper
    ) {
      stop(
        name,
        " must be one numeric value between ",
        lower,
        " and ",
        upper,
        ".",
        call. = FALSE
      )
    }
  }

  validate_number(
    single_cell_width,
    "single_cell_width",
    0.30,
    1.00
  )
  validate_number(
    single_cell_height,
    "single_cell_height",
    0.30,
    1.00
  )
  validate_number(
    multi_cell_width,
    "multi_cell_width",
    0.30,
    1.00
  )
  validate_number(
    multi_cell_height,
    "multi_cell_height",
    0.30,
    1.00
  )
  validate_number(
    figure_width,
    "figure_width",
    6.00,
    12.00
  )
  validate_number(
    figure_height,
    "figure_height",
    9.00,
    30.00
  )
  validate_number(
    trait_text_pt,
    "trait_text_pt",
    9.00,
    16.00
  )

  stop_invalid_files <- function(file_paths) {
    stop(
      "Invalid files exist: ",
      paste(basename(file_paths), collapse = ", "),
      call. = FALSE
    )
  }

  list_and_validate <- function(directory, patterns) {
    entries <- list.files(
      directory,
      all.files = TRUE,
      no.. = TRUE,
      full.names = TRUE,
      recursive = FALSE
    )

    if (length(entries) == 0L) {
      stop(
        "No result files were found in: ",
        directory,
        call. = FALSE
      )
    }

    entry_information <- file.info(entries)
    entry_names <- basename(entries)

    valid_name <- vapply(
      entry_names,
      function(entry_name) {
        any(vapply(
          patterns,
          grepl,
          x = entry_name,
          FUN.VALUE = logical(1)
        ))
      },
      FUN.VALUE = logical(1)
    )

    invalid_entries <- (
      is.na(entry_information$isdir) |
        entry_information$isdir |
        !valid_name
    )

    if (any(invalid_entries)) {
      stop_invalid_files(entries[invalid_entries])
    }

    sort(entries)
  }

  single_pattern <- paste0(
    "^PHY_kfoldcv_.+_(BM|lambda)_(FTF|TTF)[.]rds$"
  )

  all8_pattern <- paste0(
    "^PHY_kfoldcv_all8_(BM|lambda)_(FTF|TTF)[.]rds$"
  )

  cluster_pattern <- paste0(
    "^PHY_kfoldcv_high_correlation_traits_",
    "([1-6])_(BM|lambda)_(FTF|TTF)[.]rds$"
  )

  cluster_single_files <- list_and_validate(
    single_dir,
    single_pattern
  )

  all8_single_files <- list_and_validate(
    all8_single_dir,
    single_pattern
  )

  all8_files <- list_and_validate(
    all8_dir,
    all8_pattern
  )

  cluster_files <- list_and_validate(
    high_correlation_dir,
    cluster_pattern
  )

  forbidden_single_names <- grepl(
    paste0(
      "^PHY_kfoldcv_",
      "(all8|all_traits|high_correlation_traits)_"
    ),
    basename(c(cluster_single_files, all8_single_files))
  )

  if (any(forbidden_single_names)) {
    stop_invalid_files(
      c(cluster_single_files, all8_single_files)[
        forbidden_single_names
      ]
    )
  }

  expected_configurations <- c(
    "BM_FTF",
    "BM_TTF",
    "lambda_FTF",
    "lambda_TTF"
  )

  all8_traits <- c(
    "Wood_density",
    "Specific_leaf_area",
    "Seed_dry_mass",
    "Leaf_P_per_mass",
    "Stem_conduit_diameter",
    "Tree_height",
    "Root_depth",
    "Bark_thickness"
  )

  cluster_traits <- list(
    `1` = c(
      "Leaf_N_per_mass",
      "Specific_leaf_area",
      "Leaf_thickness"
    ),
    `2` = c(
      "Stem_conduit_diameter",
      "Leaf_Vcmax_per_dry_mass",
      "Stomatal_conductance",
      "Leaf_area"
    ),
    `3` = c(
      "Leaf_K_per_mass",
      "Leaf_P_per_mass"
    ),
    `4` = c(
      "Leaf_density",
      "Wood_density"
    ),
    `5` = c(
      "Bark_thickness",
      "Stem_diameter"
    ),
    `6` = c(
      "Crown_height",
      "Crown_diameter",
      "Tree_height"
    )
  )

  extract_configuration <- function(file_paths) {
    sub(
      "^.*_(BM|lambda)_(FTF|TTF)[.]rds$",
      "\\1_\\2",
      basename(file_paths)
    )
  }

  extract_single_trait <- function(file_paths) {
    names_without_extension <- tools::file_path_sans_ext(
      basename(file_paths)
    )

    sub(
      "_(BM|lambda)_(FTF|TTF)$",
      "",
      sub("^PHY_kfoldcv_", "", names_without_extension)
    )
  }

  validate_configurations <- function(file_paths, scope_label) {
    configurations <- extract_configuration(file_paths)

    if (
      length(configurations) != length(expected_configurations) ||
      !setequal(configurations, expected_configurations) ||
      anyDuplicated(configurations)
    ) {
      stop(
        scope_label,
        " must contain exactly one result for each of: ",
        paste(expected_configurations, collapse = ", "),
        call. = FALSE
      )
    }
  }

  validate_single_results <- function(
    file_paths,
    expected_traits,
    scope_label
  ) {
    observed_traits <- extract_single_trait(file_paths)

    if (!setequal(unique(observed_traits), expected_traits)) {
      stop(
        scope_label,
        " does not contain the expected trait set.",
        call. = FALSE
      )
    }

    for (trait_name in expected_traits) {
      validate_configurations(
        file_paths[observed_traits == trait_name],
        paste0(scope_label, " trait '", trait_name, "'")
      )
    }
  }

  validate_single_results(
    all8_single_files,
    all8_traits,
    "All8 single-trait directory"
  )

  validate_single_results(
    cluster_single_files,
    unname(unlist(cluster_traits, use.names = FALSE)),
    "Cluster-specific single-trait directory"
  )

  validate_configurations(
    all8_files,
    "All8 multi-trait directory"
  )

  cluster_file_names <- basename(cluster_files)
  cluster_ids <- sub(
    cluster_pattern,
    "\\1",
    cluster_file_names
  )

  for (cluster_id in names(cluster_traits)) {
    validate_configurations(
      cluster_files[cluster_ids == cluster_id],
      paste0("Cluster ", cluster_id, " multi-trait results")
    )
  }

  trait_to_cluster <- unlist(
    lapply(
      names(cluster_traits),
      function(cluster_id) {
        stats::setNames(
          rep(cluster_id, length(cluster_traits[[cluster_id]])),
          cluster_traits[[cluster_id]]
        )
      }
    )
  )

  make_file_table <- function() {
    all8_single_traits <- extract_single_trait(all8_single_files)
    cluster_single_traits <- extract_single_trait(cluster_single_files)

    cluster_single_ids <- unname(
      trait_to_cluster[cluster_single_traits]
    )

    if (anyNA(cluster_single_ids)) {
      stop(
        "At least one cluster-specific single trait could not be assigned.",
        call. = FALSE
      )
    }

    rbind(
      data.frame(
        file_path = all8_single_files,
        scope = "single",
        group_id = "all8_single",
        group_label = "8 traits",
        group_order = 0L,
        trait_name = all8_single_traits,
        configuration = extract_configuration(all8_single_files),
        stringsAsFactors = FALSE
      ),
      data.frame(
        file_path = cluster_single_files,
        scope = "single",
        group_id = paste0("cluster_", cluster_single_ids),
        group_label = as.character(
          as.roman(as.integer(cluster_single_ids))
        ),
        group_order = as.integer(cluster_single_ids),
        trait_name = cluster_single_traits,
        configuration = extract_configuration(cluster_single_files),
        stringsAsFactors = FALSE
      ),
      data.frame(
        file_path = all8_files,
        scope = "multi",
        group_id = "all8",
        group_label = "8 traits",
        group_order = 0L,
        trait_name = "8 selected traits",
        configuration = extract_configuration(all8_files),
        stringsAsFactors = FALSE
      ),
      data.frame(
        file_path = cluster_files,
        scope = "multi",
        group_id = paste0("cluster_", cluster_ids),
        group_label = as.character(
          as.roman(as.integer(cluster_ids))
        ),
        group_order = as.integer(cluster_ids),
        trait_name = paste("Cluster", cluster_ids),
        configuration = extract_configuration(cluster_files),
        stringsAsFactors = FALSE
      )
    )
  }

  file_table <- make_file_table()
  file_table$file_id <- normalizePath(
    file_table$file_path,
    winslash = "/",
    mustWork = TRUE
  )

  if (anyDuplicated(file_table$file_id)) {
    stop(
      "A result file was included more than once.",
      call. = FALSE
    )
  }

  read_fold_status <- function(file_metadata) {
    file_path <- file_metadata$file_path[[1L]]

    result <- tryCatch(
      readRDS(file_path),
      error = function(error) error
    )

    if (inherits(result, "error")) {
      stop_invalid_files(file_path)
    }

    status_table <- attr(result, "status")

    valid_structure <- (
      inherits(result, "kfold_phylopars_results") &&
        is.list(result) &&
        !is.null(names(result)) &&
        is.data.frame(status_table) &&
        all(c("fold", "status") %in% names(status_table)) &&
        nrow(status_table) == length(result)
    )

    if (!valid_structure) {
      stop_invalid_files(file_path)
    }

    fold_ids <- as.character(status_table$fold)
    status_values <- as.character(status_table$status)

    valid_status <- (
      !anyNA(fold_ids) &&
        !anyDuplicated(fold_ids) &&
        !anyNA(status_values) &&
        all(status_values %in% c("success", "failed")) &&
        all(fold_ids %in% names(result))
    )

    if (!valid_status) {
      stop_invalid_files(file_path)
    }

    result_indices <- match(fold_ids, names(result))

    has_result <- vapply(
      result[result_indices],
      function(fold_result) {
        is.list(fold_result) &&
          !is.null(fold_result$anc_recon) &&
          length(fold_result$anc_recon) > 0L &&
          !is.null(fold_result$anc_var) &&
          length(fold_result$anc_var) > 0L
      },
      FUN.VALUE = logical(1)
    )

    fold_number <- suppressWarnings(
      as.integer(gsub("[^0-9]", "", fold_ids))
    )

    if (
      anyNA(fold_number) ||
      anyDuplicated(fold_number) ||
      !setequal(fold_number, seq_len(10L))
    ) {
      stop_invalid_files(file_path)
    }

    output <- data.frame(
      file_id = file_metadata$file_id[[1L]],
      scope = file_metadata$scope[[1L]],
      group_id = file_metadata$group_id[[1L]],
      group_label = file_metadata$group_label[[1L]],
      group_order = file_metadata$group_order[[1L]],
      trait_name = file_metadata$trait_name[[1L]],
      configuration = file_metadata$configuration[[1L]],
      fold_number = fold_number,
      status = ifelse(
        status_values == "success" & has_result,
        "success",
        "failed"
      ),
      stringsAsFactors = FALSE
    )

    rm(result)
    invisible(gc(verbose = FALSE))

    output
  }

  fold_tables <- lapply(
    seq_len(nrow(file_table)),
    function(row_index) {
      read_fold_status(file_table[row_index, , drop = FALSE])
    }
  )

  fold_status <- do.call(rbind, fold_tables)
  rownames(fold_status) <- NULL

  trait_order_lookup <- c(
    stats::setNames(seq_along(all8_traits), paste0("all8_single::", all8_traits)),
    unlist(
      lapply(
        names(cluster_traits),
        function(cluster_id) {
          traits <- cluster_traits[[cluster_id]]
          stats::setNames(
            seq_along(traits),
            paste0("cluster_", cluster_id, "::", traits)
          )
        }
      )
    )
  )

  single_rows <- unique(
    file_table[
      file_table$scope == "single",
      c(
        "file_id",
        "group_id",
        "group_label",
        "group_order",
        "trait_name",
        "configuration"
      ),
      drop = FALSE
    ]
  )

  single_trait_keys <- paste0(
    single_rows$group_id,
    "::",
    single_rows$trait_name
  )

  single_rows$trait_order <- unname(
    trait_order_lookup[single_trait_keys]
  )
  single_rows$configuration_order <- match(
    single_rows$configuration,
    expected_configurations
  )

  if (
    anyNA(single_rows$trait_order) ||
    anyNA(single_rows$configuration_order)
  ) {
    stop(
      "Could not determine the single-trait row order.",
      call. = FALSE
    )
  }

  single_rows <- single_rows[
    order(
      single_rows$group_order,
      single_rows$trait_order,
      single_rows$configuration_order
    ),
    ,
    drop = FALSE
  ]

  single_rows$trait_block <- paste(
    single_rows$group_id,
    single_rows$trait_name,
    sep = "::"
  )
  block_order <- unique(single_rows$trait_block)
  single_rows$row_index <- NA_real_
  row_cursor <- 0
  previous_group <- NA_character_

  for (block_index in seq_along(block_order)) {
    block_key <- block_order[[block_index]]
    block_rows <- which(single_rows$trait_block == block_key)
    current_group <- single_rows$group_id[block_rows[[1L]]]

    if (block_index > 1L) {
      row_cursor <- row_cursor + if (
        identical(current_group, previous_group)
      ) {
        0.45
      } else {
        1.10
      }
    }

    single_rows$row_index[block_rows] <- (
      row_cursor + seq_along(block_rows)
    )
    row_cursor <- max(single_rows$row_index[block_rows])
    previous_group <- current_group
  }

  single_y_max <- max(single_rows$row_index)

  single_status <- merge(
    fold_status[fold_status$scope == "single", , drop = FALSE],
    single_rows[, c("file_id", "row_index"), drop = FALSE],
    by = "file_id",
    all.x = TRUE,
    sort = FALSE
  )

  if (anyNA(single_status$row_index)) {
    stop(
      "Could not align single-trait rows with fold statuses.",
      call. = FALSE
    )
  }

  all_trait_names <- unique(single_rows$trait_name)
  trait_colors <- stats::setNames(
    grDevices::hcl.colors(
      length(all_trait_names),
      palette = "Dark 3"
    ),
    all_trait_names
  )

  group_color_values <- c(
    all8_single = "#A6A6A6",
    all8 = "#A6A6A6",
    cluster_1 = "#D89000",
    cluster_2 = "#E64291",
    cluster_3 = "#7389C6",
    cluster_4 = "#4C9B35",
    cluster_5 = "#3579B9",
    cluster_6 = "#149977"
  )

  status_colors <- c(
    success = "#2E7D32",
    failed = "#C62828"
  )

  # All status cells use the same dimensions. A fixed coordinate ratio makes
  # them compact horizontal rectangles rather than squares.
  status_cell_size <- 0.92
  status_cell_aspect <- 0.64
  trait_text_size <- trait_text_pt / 2.845276

  wrap_trait_label <- function(trait_names, width = 18L) {
    readable_names <- gsub("_", " ", trait_names, fixed = TRUE)
    vapply(
      readable_names,
      function(trait_name) {
        paste(strwrap(trait_name, width = width), collapse = "\n")
      },
      FUN.VALUE = character(1)
    )
  }

  make_single_panel <- function() {
    trait_boxes <- do.call(
      rbind,
      lapply(
        split(
          single_rows,
          paste(single_rows$group_id, single_rows$trait_name, sep = "::")
        ),
        function(trait_rows) {
          data.frame(
            trait_name = trait_rows$trait_name[[1L]],
            ymin = min(trait_rows$row_index) - 0.45,
            ymax = max(trait_rows$row_index) + 0.45,
            y = mean(range(trait_rows$row_index)),
            stringsAsFactors = FALSE
          )
        }
      )
    )

    group_boxes <- do.call(
      rbind,
      lapply(
        split(single_rows, single_rows$group_id),
        function(group_rows) {
          data.frame(
            group_id = group_rows$group_id[[1L]],
            group_label = group_rows$group_label[[1L]],
            ymin = min(group_rows$row_index) - 0.45,
            ymax = max(group_rows$row_index) + 0.45,
            y = mean(range(group_rows$row_index)),
            stringsAsFactors = FALSE
          )
        }
      )
    )

    ggplot2::ggplot() +
      ggplot2::geom_rect(
        data = group_boxes,
        ggplot2::aes(
          xmin = -8.15,
          xmax = -7.25,
          ymin = ymin,
          ymax = ymax,
          color = group_id
        ),
        fill = "grey97",
        linewidth = 0.55,
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = group_boxes,
        ggplot2::aes(
          x = -7.70,
          y = y,
          label = group_label
        ),
        size = 2.5,
        fontface = "bold",
        angle = 90,
        lineheight = 0.95
      ) +
      ggplot2::geom_rect(
        data = trait_boxes,
        ggplot2::aes(
          xmin = -7.05,
          xmax = -0.55,
          ymin = ymin,
          ymax = ymax,
          color = trait_name
        ),
        fill = "grey99",
        linewidth = 0.50,
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = trait_boxes,
        ggplot2::aes(
          x = -3.80,
          y = y,
          label = wrap_trait_label(trait_name),
          color = trait_name
        ),
        hjust = 0.5,
        vjust = 0.5,
        size = trait_text_size,
        family = "sans",
        lineheight = 0.92,
        show.legend = FALSE
      ) +
      ggplot2::geom_tile(
        data = single_rows,
        ggplot2::aes(
          x = -1.20,
          y = row_index
        ),
        width = 2.05,
        height = single_cell_height,
        fill = NA,
        color = NA,
        linewidth = 0.40
      ) +
      ggplot2::geom_text(
        data = single_rows,
        ggplot2::aes(
          x = -1.20,
          y = row_index,
          label = ""
        ),
        size = 2.9,
        color = "grey15"
      ) +
      ggplot2::geom_tile(
        data = single_status,
        ggplot2::aes(
          x = fold_number,
          y = row_index,
          fill = status
        ),
        width = status_cell_size,
        height = status_cell_size,
        color = "grey94",
        linewidth = 0.32
      ) +
      ggplot2::geom_rect(
        data = trait_boxes,
        ggplot2::aes(
          xmin = 0.52,
          xmax = 10.48,
          ymin = ymin,
          ymax = ymax
        ),
        fill = NA,
        color = NA,
        linewidth = 0,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_color_manual(
        values = c(trait_colors, group_color_values),
        guide = "none"
      ) +
      ggplot2::scale_fill_manual(
        values = status_colors,
        guide = "none"
      ) +
      ggplot2::scale_x_continuous(
        breaks = seq_len(10L),
        labels = seq_len(10L),
        limits = c(-8.30, 10.55),
        expand = c(0, 0)
      ) +
      ggplot2::scale_y_reverse(
        limits = c(single_y_max + 0.55, 0.45),
        expand = c(0, 0)
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(
        title = NULL,
        x = "Fold",
        y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(color = "grey20"),
        axis.ticks.x = ggplot2::element_line(color = "grey50"),
        plot.title = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(8, 5, 4, 5, unit = "pt")
      )
  }

  multi_group_order <- c(
    "all8",
    paste0("cluster_", seq_len(6L))
  )

  panel_titles <- c(
    all8 = "8 selected traits",
    cluster_1 = "B(I)  Trait cluster I",
    cluster_2 = "B(II)  Trait cluster II",
    cluster_3 = "B(III)  Trait cluster III",
    cluster_4 = "B(IV)  Trait cluster IV",
    cluster_5 = "B(V)  Trait cluster V",
    cluster_6 = "B(VI)  Trait cluster VI"
  )

  make_multi_panel <- function(group_id, show_x_axis = FALSE) {
    group_rows <- unique(
      file_table[
        file_table$scope == "multi" &
          file_table$group_id == group_id,
        c("file_id", "group_id", "group_label", "configuration"),
        drop = FALSE
      ]
    )

    group_rows$configuration_order <- match(
      group_rows$configuration,
      expected_configurations
    )
    group_rows <- group_rows[
      order(group_rows$configuration_order),
      ,
      drop = FALSE
    ]

    display_traits <- if (identical(group_id, "all8")) {
      all8_traits
    } else {
      cluster_id <- sub("^cluster_", "", group_id)
      cluster_traits[[cluster_id]]
    }

    display_trait_labels <- wrap_trait_label(display_traits)
    trait_line_count <- sum(
      vapply(
        strsplit(display_trait_labels, "\n", fixed = TRUE),
        length,
        FUN.VALUE = integer(1)
      )
    )
    panel_height <- max(4, trait_line_count * 1.10)
    content_offset <- 0
    row_offset <- (panel_height - nrow(group_rows)) / 2
    group_rows$row_index <- (
      content_offset + row_offset + seq_len(nrow(group_rows))
    )

    group_status <- merge(
      fold_status[
        fold_status$scope == "multi" &
          fold_status$group_id == group_id,
        ,
        drop = FALSE
      ],
      group_rows[, c("file_id", "row_index"), drop = FALSE],
      by = "file_id",
      all.x = TRUE,
      sort = FALSE
    )

    group_traits_label <- paste(display_trait_labels, collapse = "\n")

    analysis_box <- data.frame(
      group_id = group_id,
      group_label = group_traits_label,
      ymin = content_offset + 0.55,
      ymax = content_offset + panel_height + 0.45,
      y = content_offset + (panel_height + 1) / 2,
      stringsAsFactors = FALSE
    )

    heatmap_box <- data.frame(
      ymin = content_offset + row_offset + 0.55,
      ymax = content_offset + row_offset + nrow(group_rows) + 0.45,
      stringsAsFactors = FALSE
    )

    side_label <- data.frame(
      group_id = group_id,
      label = if (identical(group_id, "all8")) {
        "8 traits"
      } else {
        cluster_number <- as.integer(sub("^cluster_", "", group_id))
        as.character(as.roman(cluster_number))
      },
      ymin = analysis_box$ymin,
      ymax = analysis_box$ymax,
      y = analysis_box$y,
      stringsAsFactors = FALSE
    )

    ggplot2::ggplot() +
      ggplot2::geom_rect(
        data = side_label,
        ggplot2::aes(
          xmin = -8.00,
          xmax = -7.15,
          ymin = ymin,
          ymax = ymax,
          color = group_id
        ),
        fill = "grey97",
        linewidth = 0.55,
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = side_label,
        ggplot2::aes(x = -7.575, y = y, label = label),
        angle = 90,
        family = "sans",
        fontface = "bold",
        size = 2.15,
        lineheight = 0.82,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_rect(
        data = analysis_box,
        ggplot2::aes(
          xmin = -6.95,
          xmax = -0.55,
          ymin = ymin,
          ymax = ymax,
          color = group_id
        ),
        fill = "grey98",
        linewidth = 0.65,
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = analysis_box,
        ggplot2::aes(
          x = -3.75,
          y = y,
          label = group_label
        ),
        hjust = 0.5,
        vjust = 0.5,
        size = trait_text_size,
        family = "sans",
        fontface = "plain",
        lineheight = 0.72
      ) +
      ggplot2::geom_tile(
        data = group_rows,
        ggplot2::aes(
          x = -0.70,
          y = row_index
        ),
        width = 1.75,
        height = multi_cell_height,
        fill = NA,
        color = NA,
        linewidth = 0.45
      ) +
      ggplot2::geom_text(
        data = group_rows,
        ggplot2::aes(
          x = -0.70,
          y = row_index,
          label = ""
        ),
        size = 3.6,
        color = "grey15"
      ) +
      ggplot2::geom_tile(
        data = group_status,
        ggplot2::aes(
          x = fold_number,
          y = row_index,
          fill = status
        ),
        width = status_cell_size,
        height = status_cell_size,
        color = "grey94",
        linewidth = 0.42
      ) +
      ggplot2::geom_rect(
        data = heatmap_box,
        ggplot2::aes(
          xmin = 0.52,
          xmax = 10.48,
          ymin = ymin,
          ymax = ymax
        ),
        fill = NA,
        color = NA,
        linewidth = 0,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_color_manual(
        values = group_color_values,
        guide = "none"
      ) +
      ggplot2::scale_fill_manual(
        values = status_colors,
        guide = "none"
      ) +
      ggplot2::scale_x_continuous(
        breaks = if (show_x_axis) seq_len(10L) else NULL,
        labels = if (show_x_axis) seq_len(10L) else NULL,
        limits = c(-8.10, 10.55),
        expand = c(0, 0)
      ) +
      ggplot2::scale_y_reverse(
        limits = c(panel_height + 0.55, 0.45),
        expand = c(0, 0)
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(
        title = NULL,
        x = if (show_x_axis) "Fold" else NULL,
        y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.x = if (show_x_axis) {
          ggplot2::element_text(color = "grey20")
        } else {
          ggplot2::element_blank()
        },
        axis.ticks.x = if (show_x_axis) {
          ggplot2::element_line(color = "grey50")
        } else {
          ggplot2::element_blank()
        },
        axis.title.x = if (show_x_axis) {
          ggplot2::element_text(color = "grey20")
        } else {
          ggplot2::element_blank()
        },
        plot.title = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(8, 5, 4, 5, unit = "pt")
      )
  }

  make_multi_header_panel <- function() {
    ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0.5,
        label = "B  Multivariate analyses",
        hjust = 0,
        vjust = 0.5,
        fontface = "bold",
        size = 5.2
      ) +
      ggplot2::coord_cartesian(
        xlim = c(0, 1),
        ylim = c(0, 1),
        clip = "off"
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0, 5, 0, 5, unit = "pt")
      )
  }

  make_column_header_panel <- function(label) {
    ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0.5,
        label = label,
        hjust = 0,
        vjust = 0.5,
        family = "sans",
        fontface = "bold",
        size = 14 / 2.845276
      ) +
      ggplot2::coord_cartesian(
        xlim = c(0, 1),
        ylim = c(0, 1),
        clip = "off"
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0, 5, 0, 5, unit = "pt")
      )
  }

  make_legend_panel <- function() {
    ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0.02,
        y = 13.6,
        label = "Status",
        hjust = 0,
        fontface = "bold",
        size = 5.2
      ) +
      ggplot2::annotate(
        "rect",
        xmin = 0.02,
        xmax = 0.08,
        ymin = 12.75,
        ymax = 13.15,
        fill = status_colors[["success"]],
        color = NA
      ) +
      ggplot2::annotate(
        "text",
        x = 0.10,
        y = 12.95,
        label = "Success",
        hjust = 0,
        size = 4.6
      ) +
      ggplot2::annotate(
        "rect",
        xmin = 0.02,
        xmax = 0.08,
        ymin = 12.05,
        ymax = 12.45,
        fill = status_colors[["failed"]],
        color = NA
      ) +
      ggplot2::annotate(
        "text",
        x = 0.10,
        y = 12.25,
        label = "Failed",
        hjust = 0,
        size = 4.6
      ) +
      ggplot2::annotate(
        "segment",
        x = 0.02,
        xend = 0.96,
        y = 11.55,
        yend = 11.55,
        color = "grey82",
        linewidth = 0.45
      ) +
      ggplot2::annotate(
        "text",
        x = 0.02,
        y = 11.05,
        label = "Row order within each four-row block",
        hjust = 0,
        fontface = "bold",
        size = 5.2
      ) +
      ggplot2::annotate(
        "text",
        x = 0.02,
        y = 10.25,
        label = paste0(
          "1  BM-FTF\n",
          "2  BM-TTF\n",
          "3  lambda-FTF\n",
          "4  lambda-TTF"
        ),
        hjust = 0,
        vjust = 1,
        fontface = "plain",
        lineheight = 1.10,
        size = 4.2
      ) +
      ggplot2::annotate(
        "text",
        x = 0.18,
        y = 10.25,
        label = "",
        hjust = 0,
        vjust = 1,
        lineheight = 1.15,
        size = 4.4
      ) +
      ggplot2::annotate(
        "text",
        x = 0.02,
        y = 7.45,
        label = "Configuration codes",
        hjust = 0,
        vjust = 1,
        fontface = "bold",
        size = 4.7
      ) +
      ggplot2::annotate(
        "text",
        x = 0.02,
        y = 6.75,
        label = paste0(
          "FTF: pheno_error = FALSE\n",
          "TTF: pheno_error = TRUE\n",
          "Both: phylo_correlated = TRUE\n",
          "      pheno_correlated = FALSE"
        ),
        hjust = 0,
        vjust = 1,
        lineheight = 1.15,
        size = 4.4
      ) +
      ggplot2::annotate(
        "segment",
        x = 0.02,
        xend = 0.96,
        y = 3.95,
        yend = 3.95,
        color = "grey82",
        linewidth = 0.45
      ) +
      ggplot2::annotate(
        "text",
        x = 0.02,
        y = 3.45,
        label = "Evolutionary models",
        hjust = 0,
        fontface = "bold",
        size = 5.2
      ) +
      ggplot2::annotate(
        "text",
        x = 0.02,
        y = 2.55,
        label = "BM",
        hjust = 0,
        fontface = "bold",
        size = 4.7
      ) +
      ggplot2::annotate(
        "text",
        x = 0.26,
        y = 2.55,
        label = "Brownian motion",
        hjust = 0,
        size = 4.4
      ) +
      ggplot2::annotate(
        "text",
        x = 0.02,
        y = 1.65,
        label = "lambda",
        hjust = 0,
        fontface = "bold",
        size = 4.7
      ) +
      ggplot2::annotate(
        "text",
        x = 0.26,
        y = 1.65,
        label = "Pagel's lambda",
        hjust = 0,
        size = 4.4
      ) +
      ggplot2::coord_cartesian(
        xlim = c(0, 1),
        ylim = c(0.90, 13.90),
        clip = "off"
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(5, 5, 5, 5, unit = "pt")
      )
  }

  make_compact_legend_panel <- function() {
    heading_size <- 3.65
    body_size <- 3.20

    ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0.02,
        y = 2.75,
        label = "Status",
        hjust = 0,
        fontface = "bold",
        size = heading_size
      ) +
      ggplot2::annotate(
        "rect",
        xmin = 0.02,
        xmax = 0.045,
        ymin = 2.05,
        ymax = 2.30,
        fill = status_colors[["success"]],
        color = NA
      ) +
      ggplot2::annotate(
        "text",
        x = 0.055,
        y = 2.175,
        label = "Success",
        hjust = 0,
        size = body_size
      ) +
      ggplot2::annotate(
        "rect",
        xmin = 0.02,
        xmax = 0.045,
        ymin = 1.55,
        ymax = 1.80,
        fill = status_colors[["failed"]],
        color = NA
      ) +
      ggplot2::annotate(
        "text",
        x = 0.055,
        y = 1.675,
        label = "Failed",
        hjust = 0,
        size = body_size
      ) +
      ggplot2::annotate(
        "text",
        x = 0.19,
        y = 2.75,
        label = "Row order (top to bottom)",
        hjust = 0,
        fontface = "bold",
        size = heading_size
      ) +
      ggplot2::annotate(
        "text",
        x = 0.19,
        y = 2.30,
        label = paste0(
          "1  BM-FTF\n",
          "2  BM-TTF\n",
          "3  lambda-FTF\n",
          "4  lambda-TTF"
        ),
        hjust = 0,
        vjust = 1,
        lineheight = 1.10,
        size = body_size
      ) +
      ggplot2::annotate(
        "text",
        x = 0.49,
        y = 2.75,
        label = "Configuration codes",
        hjust = 0,
        fontface = "bold",
        size = heading_size
      ) +
      ggplot2::annotate(
        "text",
        x = 0.49,
        y = 2.30,
        label = paste0(
          "FTF: pheno_error = FALSE\n",
          "TTF: pheno_error = TRUE\n",
          "Both: phylo_correlated = TRUE\n",
          "       pheno_correlated = FALSE"
        ),
        hjust = 0,
        vjust = 1,
        lineheight = 1.10,
        size = body_size
      ) +
      ggplot2::annotate(
        "text",
        x = 0.80,
        y = 2.75,
        label = "Evolutionary models",
        hjust = 0,
        fontface = "bold",
        size = heading_size
      ) +
      ggplot2::annotate(
        "text",
        x = 0.80,
        y = 2.30,
        label = "BM      Brownian motion\nlambda  Pagel's lambda",
        hjust = 0,
        vjust = 1,
        lineheight = 1.25,
        size = body_size
      ) +
      ggplot2::coord_cartesian(
        xlim = c(0, 1),
        ylim = c(0.75, 3.05),
        clip = "off"
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(5, 5, 4, 5, unit = "pt")
      )
  }

  make_b_legend_panel <- function() {
    heading_size <- 3.65
    body_size <- 3.20

    ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = 0.02, y = 12.70,
        label = "Status", hjust = 0,
        family = "sans", fontface = "bold", size = heading_size
      ) +
      ggplot2::annotate(
        "rect", xmin = 0.02, xmax = 0.075,
        ymin = 11.85, ymax = 12.15,
        fill = status_colors[["success"]], color = NA
      ) +
      ggplot2::annotate(
        "text", x = 0.095, y = 12.00,
        label = "Success", hjust = 0,
        family = "sans", size = body_size
      ) +
      ggplot2::annotate(
        "rect", xmin = 0.02, xmax = 0.075,
        ymin = 11.25, ymax = 11.55,
        fill = status_colors[["failed"]], color = NA
      ) +
      ggplot2::annotate(
        "text", x = 0.095, y = 11.40,
        label = "Failed", hjust = 0,
        family = "sans", size = body_size
      ) +
      ggplot2::annotate(
        "text", x = 0.02, y = 10.55,
        label = "Row order (top to bottom)", hjust = 0,
        family = "sans", fontface = "bold", size = heading_size
      ) +
      ggplot2::annotate(
        "text", x = 0.02, y = 10.05,
        label = paste0(
          "1  BM-FTF\n2  BM-TTF\n",
          "3  lambda-FTF\n4  lambda-TTF"
        ),
        hjust = 0, vjust = 1, lineheight = 1.00,
        family = "sans", size = body_size
      ) +
      ggplot2::annotate(
        "text", x = 0.02, y = 6.55,
        label = "Configuration codes", hjust = 0,
        family = "sans", fontface = "bold", size = heading_size
      ) +
      ggplot2::annotate(
        "text", x = 0.02, y = 6.05,
        label = paste0(
          "FTF: pheno_error = FALSE\n",
          "TTF: pheno_error = TRUE\n",
          "Both: phylo_correlated = TRUE\n",
          "       pheno_correlated = FALSE"
        ),
        hjust = 0, vjust = 1, lineheight = 1.00,
        family = "sans", size = body_size
      ) +
      ggplot2::annotate(
        "text", x = 0.02, y = 2.55,
        label = "Evolutionary models", hjust = 0,
        family = "sans", fontface = "bold", size = heading_size
      ) +
      ggplot2::annotate(
        "text", x = 0.02, y = 2.05,
        label = "BM      Brownian motion\nlambda  Pagel's lambda",
        hjust = 0, vjust = 1, lineheight = 1.25,
        family = "sans", size = body_size
      ) +
      ggplot2::coord_cartesian(
        xlim = c(0, 1), ylim = c(0.35, 13.10), clip = "off"
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(5, 5, 4, 5, unit = "pt")
      )
  }

  make_b_legend_panel_compact <- function() {
    heading_size <- 3.45
    body_size <- 3.20

    ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = 0.02, y = 5.85,
        label = "Status", hjust = 0,
        family = "sans", fontface = "bold", size = heading_size
      ) +
      ggplot2::annotate(
        "rect", xmin = 0.02, xmax = 0.07,
        ymin = 5.15, ymax = 5.42,
        fill = status_colors[["success"]], color = NA
      ) +
      ggplot2::annotate(
        "text", x = 0.085, y = 5.285,
        label = "Success", hjust = 0,
        family = "sans", size = body_size
      ) +
      ggplot2::annotate(
        "rect", xmin = 0.02, xmax = 0.07,
        ymin = 4.65, ymax = 4.92,
        fill = status_colors[["failed"]], color = NA
      ) +
      ggplot2::annotate(
        "text", x = 0.085, y = 4.785,
        label = "Failed", hjust = 0,
        family = "sans", size = body_size
      ) +
      ggplot2::annotate(
        "text", x = 0.02, y = 3.95,
        label = "Row order (top to bottom)", hjust = 0,
        family = "sans", fontface = "bold", size = heading_size
      ) +
      ggplot2::annotate(
        "text", x = 0.02, y = 3.55,
        label = paste0(
          "1  BM-FTF\n2  BM-TTF\n",
          "3  lambda-FTF\n4  lambda-TTF"
        ),
        hjust = 0, vjust = 1, lineheight = 1.00,
        family = "sans", size = body_size
      ) +
      ggplot2::annotate(
        "text", x = 0.52, y = 5.85,
        label = "Configuration codes", hjust = 0,
        family = "sans", fontface = "bold", size = heading_size
      ) +
      ggplot2::annotate(
        "text", x = 0.52, y = 5.45,
        label = paste0(
          "FTF: pheno_error = FALSE\n",
          "TTF: pheno_error = TRUE\n",
          "Both: phylo_correlated = TRUE\n",
          "       pheno_correlated = FALSE"
        ),
        hjust = 0, vjust = 1, lineheight = 1.00,
        family = "sans", size = body_size
      ) +
      ggplot2::annotate(
        "text", x = 0.52, y = 2.60,
        label = "Evolutionary models", hjust = 0,
        family = "sans", fontface = "bold", size = heading_size
      ) +
      ggplot2::annotate(
        "text", x = 0.52, y = 2.20,
        label = "BM      Brownian motion\nlambda  Pagel's lambda",
        hjust = 0, vjust = 1, lineheight = 1.15,
        family = "sans", size = body_size
      ) +
      ggplot2::coord_cartesian(
        xlim = c(0, 1), ylim = c(0.80, 6.15), clip = "off"
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(3, 5, 3, 5, unit = "pt")
      )
  }

  panel_single <- make_single_panel()

  header_left <- make_column_header_panel("A  Univariate analyses")
  header_right <- make_column_header_panel("B  Multivariate analyses")

  right_panels <- lapply(
    seq_along(multi_group_order),
    function(index) {
      make_multi_panel(
        multi_group_order[[index]],
        show_x_axis = index == length(multi_group_order)
      )
    }
  )

  legend_panel <- make_b_legend_panel_compact()

  multi_label_height <- function(traits) {
    labels <- wrap_trait_label(traits)
    line_count <- sum(
      vapply(
        strsplit(labels, "\n", fixed = TRUE),
        length,
        FUN.VALUE = integer(1)
      )
    )
    max(4, line_count * 1.10)
  }

  right_panel_heights <- c(
    multi_label_height(all8_traits),
    vapply(cluster_traits, multi_label_height, FUN.VALUE = numeric(1))
  )
  right_panel_heights <- pmax(4.5, right_panel_heights + 0.80)
  legend_height <- 12
  right_content_height <- (
    sum(right_panel_heights) +
      legend_height
  )
  main_content_height <- single_y_max * status_cell_aspect + 7
  spacer_height <- max(
    1,
    main_content_height - right_content_height
  )

  right_column <- patchwork::wrap_plots(
    c(
      right_panels,
      list(legend_panel, patchwork::plot_spacer())
    ),
    ncol = 1,
    heights = c(
      right_panel_heights,
      legend_height,
      spacer_height
    )
  )

  header_row <- patchwork::wrap_plots(
    header_left,
    header_right,
    ncol = 2,
    widths = c(1, 1)
  )

  content_row <- patchwork::wrap_plots(
    panel_single,
    right_column,
    ncol = 2,
    widths = c(1, 1)
  )

  main_panels <- patchwork::wrap_plots(
    header_row,
    content_row,
    ncol = 1,
    heights = c(3, main_content_height)
  )

  status_plot <- main_panels +
    patchwork::plot_annotation(
      title = "Fold Status",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          size = 19,
          hjust = 0,
          margin = ggplot2::margin(b = 5, unit = "pt")
        ),
        plot.margin = ggplot2::margin(4, 0, 0, 0, unit = "pt")
      )
    )

  if (!is.null(output_file)) {
    if (
      length(output_file) != 1L ||
      is.na(output_file) ||
      !nzchar(trimws(output_file))
    ) {
      stop(
        "Exactly one output file must be provided.",
        call. = FALSE
      )
    }

    output_directory <- dirname(output_file)

    if (!dir.exists(output_directory)) {
      stop(
        "Output directory does not exist.",
        call. = FALSE
      )
    }

    output_extension <- tolower(tools::file_ext(output_file))

    if (!output_extension %in% c("png", "pdf")) {
      stop(
        "Output format must be PNG or PDF.",
        call. = FALSE
      )
    }

    save_arguments <- list(
      filename = output_file,
      plot = status_plot,
      width = figure_width,
      height = figure_height,
      units = "in",
      dpi = 300,
      bg = "white"
    )

    if (
      output_extension == "png" &&
      requireNamespace("ragg", quietly = TRUE)
    ) {
      save_arguments$device <- ragg::agg_png
    }

    do.call(ggplot2::ggsave, save_arguments)
  }

  print(status_plot)

  invisible(status_plot)
}

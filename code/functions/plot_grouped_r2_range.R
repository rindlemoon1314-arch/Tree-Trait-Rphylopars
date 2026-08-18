plot_grouped_r2_range <- function(
    single_dir = "",
    all8_single_dir = "",
    all8_dir = "",
    high_correlation_dir = "",
    folds = c(1, 3, 5, 6, 8, 10),
    negative_compression = 0.05,
    output_file = NULL,
    data_file = NULL,
    base_font_size = 9,
    title_font_size = 10.5,
    output_width = 8.25,
    minimum_height = 10.5,
    row_height = 0.27,
    font_family = "sans"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required.",
      call. = FALSE
    )
  }

  if (!exists("calculate_phy_metrics", mode = "function")) {
    stop(
      paste0(
        "calculate_phy_metrics() is not available. ",
        "Source calculate_phy_metrics.R first."
      ),
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

  if (
    !is.numeric(folds) ||
    length(folds) == 0L ||
    anyNA(folds) ||
    any(!is.finite(folds)) ||
    any(folds < 1) ||
    any(folds != as.integer(folds)) ||
    anyDuplicated(folds)
  ) {
    stop(
      "folds must contain unique positive integers.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(negative_compression) ||
    length(negative_compression) != 1L ||
    is.na(negative_compression) ||
    !is.finite(negative_compression) ||
    negative_compression <= 0 ||
    negative_compression > 1
  ) {
    stop(
      "negative_compression must be greater than 0 and no greater than 1.",
      call. = FALSE
    )
  }

  configurations <- c(
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
    `Cluster 1` = c(
      "Leaf_N_per_mass",
      "Specific_leaf_area",
      "Leaf_thickness"
    ),
    `Cluster 2` = c(
      "Stem_conduit_diameter",
      "Leaf_Vcmax_per_dry_mass",
      "Stomatal_conductance",
      "Leaf_area"
    ),
    `Cluster 3` = c(
      "Leaf_K_per_mass",
      "Leaf_P_per_mass"
    ),
    `Cluster 4` = c(
      "Leaf_density",
      "Wood_density"
    ),
    `Cluster 5` = c(
      "Bark_thickness",
      "Stem_diameter"
    ),
    `Cluster 6` = c(
      "Crown_height",
      "Crown_diameter",
      "Tree_height"
    )
  )

  group_traits <- c(
    list(All8 = all8_traits),
    cluster_traits
  )

  extract_configuration <- function(file_paths) {
    sub(
      "^.*_(BM|lambda)_(FTF|TTF)[.]rds$",
      "\\1_\\2",
      basename(file_paths)
    )
  }

  extract_single_trait <- function(file_paths) {
    file_names <- tools::file_path_sans_ext(
      basename(file_paths)
    )

    sub(
      "_(BM|lambda)_(FTF|TTF)$",
      "",
      sub("^PHY_kfoldcv_", "", file_names)
    )
  }

  validate_configurations <- function(file_paths, label) {
    observed <- extract_configuration(file_paths)

    if (
      length(observed) != length(configurations) ||
      !setequal(observed, configurations) ||
      anyDuplicated(observed)
    ) {
      stop(
        label,
        " must contain exactly one result for each of: ",
        paste(configurations, collapse = ", "),
        call. = FALSE
      )
    }
  }

  validate_single_files <- function(
    file_paths,
    expected_traits,
    label
  ) {
    observed_traits <- extract_single_trait(file_paths)

    if (!setequal(unique(observed_traits), expected_traits)) {
      stop(
        label,
        " does not contain the expected trait set.",
        call. = FALSE
      )
    }

    for (trait_name in expected_traits) {
      validate_configurations(
        file_paths[observed_traits == trait_name],
        paste0(label, " trait '", trait_name, "'")
      )
    }
  }

  single_pattern <- paste0(
    "^PHY_kfoldcv_.+_(BM|lambda)_(FTF|TTF)[.]rds$"
  )

  cluster_single_files <- list.files(
    single_dir,
    pattern = single_pattern,
    full.names = TRUE
  )

  all8_single_files <- list.files(
    all8_single_dir,
    pattern = single_pattern,
    full.names = TRUE
  )

  all8_multi_files <- list.files(
    all8_dir,
    pattern = paste0(
      "^PHY_kfoldcv_all8_",
      "(BM|lambda)_(FTF|TTF)[.]rds$"
    ),
    full.names = TRUE
  )

  cluster_multi_files <- list.files(
    high_correlation_dir,
    pattern = paste0(
      "^PHY_kfoldcv_high_correlation_traits_",
      "([1-6])_(BM|lambda)_(FTF|TTF)[.]rds$"
    ),
    full.names = TRUE
  )

  validate_single_files(
    all8_single_files,
    all8_traits,
    "All8 single-trait results"
  )

  validate_single_files(
    cluster_single_files,
    unname(unlist(cluster_traits, use.names = FALSE)),
    "Cluster-specific single-trait results"
  )

  validate_configurations(
    all8_multi_files,
    "All8 multi-trait results"
  )

  cluster_multi_names <- basename(cluster_multi_files)
  cluster_ids <- sub(
    paste0(
      "^PHY_kfoldcv_high_correlation_traits_",
      "([1-6])_(BM|lambda)_(FTF|TTF)[.]rds$"
    ),
    "\\1",
    cluster_multi_names
  )

  for (cluster_id in seq_len(6L)) {
    validate_configurations(
      cluster_multi_files[cluster_ids == as.character(cluster_id)],
      paste("Cluster", cluster_id, "multi-trait results")
    )
  }

  cluster_single_trait_names <- extract_single_trait(
    cluster_single_files
  )

  group_files <- vector("list", length(group_traits))
  names(group_files) <- names(group_traits)

  group_files[["All8"]] <- c(
    all8_single_files,
    all8_multi_files
  )

  for (cluster_id in seq_len(6L)) {
    group_label <- paste("Cluster", cluster_id)
    traits <- cluster_traits[[group_label]]

    group_files[[group_label]] <- c(
      cluster_single_files[
        cluster_single_trait_names %in% traits
      ],
      cluster_multi_files[
        cluster_ids == as.character(cluster_id)
      ]
    )
  }

  result_tables <- vector("list", length(group_files))

  for (group_index in seq_along(group_files)) {
    group_label <- names(group_files)[[group_index]]
    current_files <- group_files[[group_index]]

    names(current_files) <- tools::file_path_sans_ext(
      basename(current_files)
    )

    if (anyDuplicated(names(current_files))) {
      stop(
        "Duplicate result names were found in ",
        group_label,
        ".",
        call. = FALSE
      )
    }

    message(
      sprintf(
        "[%d/%d] Calculating R2: %s",
        group_index,
        length(group_files),
        group_label
      )
    )

    current_metrics <- calculate_phy_metrics(
      phy_results = current_files,
      traits = group_traits[[group_label]],
      intervals = 95,
      use_transformed_truth = TRUE,
      folds = folds,
      save_dir = NULL,
      save_name = paste0(
        "r2_",
        gsub(" ", "_", tolower(group_label), fixed = TRUE)
      )
    )

    current_table <- current_metrics$trait_metrics

    is_multi <- if (identical(group_label, "All8")) {
      grepl("^PHY_kfoldcv_all8_", current_table$result_name)
    } else {
      cluster_number <- sub("^Cluster ", "", group_label)
      grepl(
        paste0(
          "^PHY_kfoldcv_high_correlation_traits_",
          cluster_number,
          "_"
        ),
        current_table$result_name
      )
    }

    current_table$group <- group_label
    current_table$analysis <- ifelse(
      is_multi,
      "Multi-trait",
      "Single trait"
    )
    current_table$configuration_label <- paste(
      current_table$model,
      current_table$configuration,
      sep = "_"
    )

    expected_rows_per_trait <- 8L
    observed_counts <- table(current_table$trait)

    if (
      !setequal(
        unique(current_table$trait),
        group_traits[[group_label]]
      ) ||
      any(observed_counts != expected_rows_per_trait)
    ) {
      stop(
        "The R2 result grid is incomplete for ",
        group_label,
        ".",
        call. = FALSE
      )
    }

    result_tables[[group_index]] <- current_table[
      ,
      c(
        "group",
        "trait",
        "analysis",
        "model",
        "configuration",
        "configuration_label",
        "folds_used",
        "r2"
      ),
      drop = FALSE
    ]

    rm(current_metrics, current_table)
    invisible(gc(verbose = FALSE))

    message(
      sprintf(
        "[%d/%d] R2 completed: %s",
        group_index,
        length(group_files),
        group_label
      )
    )
  }

  plot_data <- do.call(rbind, result_tables)
  rownames(plot_data) <- NULL

  group_order <- names(group_traits)

  row_tables <- vector("list", length(group_order))
  cursor <- 1

  for (group_index in seq_along(group_order)) {
    group_label <- group_order[[group_index]]
    traits <- group_traits[[group_label]]

    row_tables[[group_index]] <- data.frame(
      group = group_label,
      trait = traits,
      group_order = group_index,
      trait_order = seq_along(traits),
      row_position = cursor + seq_along(traits) - 1,
      stringsAsFactors = FALSE
    )

    cursor <- cursor + length(traits) + 1
  }

  row_table <- do.call(rbind, row_tables)
  rownames(row_table) <- NULL

  plot_data <- merge(
    plot_data,
    row_table[, c("group", "trait", "row_position")],
    by = c("group", "trait"),
    all.x = TRUE,
    sort = FALSE
  )

  if (anyNA(plot_data$row_position)) {
    stop(
      "Could not align R2 results with plot rows.",
      call. = FALSE
    )
  }

  configuration_offsets <- c(
    BM_FTF = -0.18,
    BM_TTF = -0.06,
    lambda_FTF = 0.06,
    lambda_TTF = 0.18
  )

  analysis_offsets <- c(
    `Single trait` = -0.018,
    `Multi-trait` = 0.018
  )

  plot_data$point_position <- (
    plot_data$row_position +
      unname(configuration_offsets[plot_data$configuration_label]) +
      unname(analysis_offsets[plot_data$analysis])
  )

  if (anyNA(plot_data$point_position)) {
    stop(
      "Could not determine point positions.",
      call. = FALSE
    )
  }

  range_keys <- interaction(
    plot_data$group,
    plot_data$trait,
    drop = TRUE,
    lex.order = TRUE
  )

  range_rows <- lapply(
    split(seq_len(nrow(plot_data)), range_keys),
    function(indices) {
      current <- plot_data[indices, , drop = FALSE]

      data.frame(
        group = current$group[[1L]],
        trait = current$trait[[1L]],
        row_position = current$row_position[[1L]],
        minimum_r2 = min(current$r2, na.rm = TRUE),
        maximum_r2 = max(current$r2, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  )

  range_data <- do.call(rbind, range_rows)
  rownames(range_data) <- NULL

  compress_negative_values <- function(values) {
    ifelse(
      values < 0,
      values * negative_compression,
      values
    )
  }

  plot_data$plot_r2 <- compress_negative_values(plot_data$r2)
  range_data$minimum_plot_r2 <- compress_negative_values(
    range_data$minimum_r2
  )
  range_data$maximum_plot_r2 <- compress_negative_values(
    range_data$maximum_r2
  )

  negative_data <- plot_data[
    is.finite(plot_data$r2) & plot_data$r2 < 0,
    ,
    drop = FALSE
  ]

  negative_data$r2_label <- formatC(
    negative_data$r2,
    format = "f",
    digits = 2
  )
  negative_data$label_hjust <- ifelse(
    negative_data$r2 > -0.5,
    1.18,
    -0.18
  )
  negative_data$label_position <- negative_data$point_position

  duplicate_label_groups <- interaction(
    negative_data$group,
    negative_data$trait,
    drop = TRUE
  )

  for (indices in split(seq_len(nrow(negative_data)), duplicate_label_groups)) {
    if (length(indices) > 1L) {
      center <- unique(negative_data$row_position[indices])[[1L]]
      negative_data$label_position[indices] <- seq(
        center - 0.20,
        center + 0.20,
        length.out = length(indices)
      )
    }
  }

  group_colors <- c(
    All8 = "#4D4D4D",
    `Cluster 1` = "#D89000",
    `Cluster 2` = "#E64291",
    `Cluster 3` = "#7389C6",
    `Cluster 4` = "#4C9B35",
    `Cluster 5` = "#3579B9",
    `Cluster 6` = "#149977"
  )

  configuration_colors <- c(
    BM_FTF = "#0072B2",
    BM_TTF = "#E69F00",
    lambda_FTF = "#56B4E9",
    lambda_TTF = "#D55E00"
  )

  range_data$group_color <- unname(
    group_colors[range_data$group]
  )

  group_label_data <- do.call(
    rbind,
    lapply(
      group_order,
      function(group_label) {
        positions <- row_table$row_position[
          row_table$group == group_label
        ]

        display_group <- if (identical(group_label, "All8")) {
          "8 selected\ntraits"
        } else {
          group_label
        }

        data.frame(
          group = group_label,
          display_group = display_group,
          row_position = mean(range(positions)),
          group_color = unname(group_colors[[group_label]]),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  separator_positions <- vapply(
    group_order[-length(group_order)],
    function(group_label) {
      max(row_table$row_position[row_table$group == group_label]) + 1
    },
    numeric(1)
  )

  maximum_r2 <- max(plot_data$r2, na.rm = TRUE)

  x_minimum <- min(
    0,
    min(plot_data$plot_r2, na.rm = TRUE) * 1.08
  )
  x_maximum <- max(0.1, ceiling(maximum_r2 * 10) / 10)

  positive_breaks <- seq(
    0,
    x_maximum + sqrt(.Machine$double.eps),
    by = 0.1
  )

  axis_breaks <- positive_breaks

  plot_data$configuration_label <- factor(
    plot_data$configuration_label,
    levels = configurations
  )
  plot_data$analysis <- factor(
    plot_data$analysis,
    levels = c("Single trait", "Multi-trait")
  )

  r2_plot <- ggplot2::ggplot() +
    ggplot2::geom_hline(
      yintercept = separator_positions,
      color = "grey68",
      linewidth = 0.55,
      linetype = "dashed"
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      color = "grey42",
      linewidth = 0.65
    ) +
    ggplot2::geom_segment(
      data = range_data,
      ggplot2::aes(
        x = minimum_plot_r2,
        xend = maximum_plot_r2,
        y = row_position,
        yend = row_position
      ),
      color = range_data$group_color,
      linewidth = 1.05,
      alpha = 0.72
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(
        x = plot_r2,
        y = point_position,
        color = configuration_label,
        shape = analysis
      ),
      size = 2.8,
      stroke = 0.95
    ) +
    ggplot2::geom_segment(
      data = negative_data,
      ggplot2::aes(
        x = plot_r2,
        xend = plot_r2,
        y = point_position,
        yend = label_position
      ),
      color = "grey55",
      linewidth = 0.28,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = negative_data,
      ggplot2::aes(
        x = plot_r2,
        y = label_position,
        label = r2_label,
        color = configuration_label
      ),
      hjust = negative_data$label_hjust,
      vjust = 0.5,
      size = (base_font_size * 0.88) / 2.845276,
      fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = group_label_data,
      ggplot2::aes(
        x = Inf,
        y = row_position,
        label = display_group
      ),
      color = group_label_data$group_color,
      hjust = -0.12,
      fontface = "bold",
      size = (base_font_size * 0.98) / 2.845276
    ) +
    ggplot2::scale_color_manual(
      values = configuration_colors,
      breaks = configurations,
      labels = c(
        "BM · FTF",
        "BM · TTF",
        "lambda · FTF",
        "lambda · TTF"
      ),
      name = "Configuration"
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        `Single trait` = 1,
        `Multi-trait` = 16
      ),
      labels = c(
        `Single trait` = "Univariate",
        `Multi-trait` = "Multivariate"
      ),
      name = "Analysis"
    ) +
    ggplot2::scale_x_continuous(
      breaks = axis_breaks,
      labels = function(values) {
        formatC(values, format = "f", digits = 1)
      },
      limits = c(x_minimum, x_maximum),
      expand = ggplot2::expansion(mult = c(0.015, 0.025))
    ) +
    ggplot2::scale_y_reverse(
      breaks = row_table$row_position,
      labels = gsub(
        "_",
        " ",
        row_table$trait,
        fixed = TRUE
      ),
      limits = c(max(row_table$row_position) + 0.6, 0.4),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = "Trait-level R²",
      subtitle = paste0(
        "R² calculated across folds ",
        paste(folds, collapse = ", "),
        "\nNegative R² positions are compressed and labelled directly"
      ),
      x = "R²",
      y = NULL
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        order = 1,
        override.aes = list(shape = 16)
      ),
      shape = ggplot2::guide_legend(
        order = 2,
        override.aes = list(color = "grey25")
      )
    ) +
    ggplot2::theme_minimal(
      base_size = base_font_size,
      base_family = font_family
    ) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(
        color = "grey88",
        linewidth = 0.45
      ),
      axis.text.y = ggplot2::element_text(
        color = "grey18",
        size = base_font_size * 0.95
      ),
      axis.text.x = ggplot2::element_text(
        color = "grey18",
        size = base_font_size * 0.90
      ),
      axis.title.x = ggplot2::element_text(
        color = "grey15",
        size = base_font_size,
        margin = ggplot2::margin(t = 8)
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = title_font_size
      ),
      plot.subtitle = ggplot2::element_text(
        color = "grey38",
        size = base_font_size * 0.95,
        margin = ggplot2::margin(b = 10)
      ),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(
        face = "bold",
        size = base_font_size * 0.95
      ),
      legend.text = ggplot2::element_text(
        size = base_font_size * 0.90
      ),
      plot.margin = ggplot2::margin(
        10,
        95,
        10,
        10,
        unit = "pt"
      )
    )

  if (!is.null(data_file)) {
    if (
      !is.character(data_file) ||
      length(data_file) != 1L ||
      is.na(data_file) ||
      !nzchar(trimws(data_file))
    ) {
      stop(
        "data_file must be NULL or one non-empty path.",
        call. = FALSE
      )
    }

    dir.create(
      dirname(data_file),
      recursive = TRUE,
      showWarnings = FALSE
    )

    utils::write.csv(
      plot_data[
        order(
          match(plot_data$group, group_order),
          plot_data$row_position,
          match(plot_data$configuration_label, configurations),
          match(plot_data$analysis, c("Single trait", "Multi-trait"))
        ),
        c(
          "group",
          "trait",
          "analysis",
          "configuration_label",
          "folds_used",
          "r2"
        )
      ],
      data_file,
      row.names = FALSE,
      na = ""
    )
  }

  if (!is.null(output_file)) {
    if (
      !is.character(output_file) ||
      length(output_file) != 1L ||
      is.na(output_file) ||
      !nzchar(trimws(output_file))
    ) {
      stop(
        "output_file must be NULL or one non-empty path.",
        call. = FALSE
      )
    }

    extension <- tolower(tools::file_ext(output_file))

    if (!extension %in% c("png", "pdf")) {
      stop(
        "The output format must be PNG or PDF.",
        call. = FALSE
      )
    }

    dir.create(
      dirname(output_file),
      recursive = TRUE,
      showWarnings = FALSE
    )

    save_arguments <- list(
      filename = output_file,
      plot = r2_plot,
      width = output_width,
      height = max(
        minimum_height,
        row_height * nrow(row_table) + 2.8
      ),
      units = "in",
      dpi = 300,
      bg = "white"
    )

    if (
      extension == "png" &&
      requireNamespace("ragg", quietly = TRUE)
    ) {
      save_arguments$device <- ragg::agg_png
    }

    do.call(ggplot2::ggsave, save_arguments)
  }

  print(r2_plot)

  invisible(
    list(
      plot = r2_plot,
      r2_data = plot_data,
      range_data = range_data,
      row_table = row_table
    )
  )
}

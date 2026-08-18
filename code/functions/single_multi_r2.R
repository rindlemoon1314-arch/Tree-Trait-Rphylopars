plot_single_multi_r2_comparison <- function(
    r2_data,
    output_file = NULL,
    heatmap_file = NULL,
    combined_file = NULL,
    detailed_table_file = NULL,
    summary_table_file = NULL,
    negative_compression = 0.035,
    delta_fill_limit = 0.05,
    width = 6.5,
    height = 9.0,
    combined_width = 8.25,
    combined_height = 11.0,
    dpi = 300,
    base_font_size = 7.2,
    title_font_size = 8.5,
    panel_title_font_size = 13.5,
    trait_spacing_multiplier = 1.00,
    font_family = "Arial"
) {
  required_packages <- c("ggplot2", "patchwork", "scales")
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

  # ggplot2 geom text sizes are expressed in millimetres.  These helpers
  # keep labels proportional to the manuscript-oriented point size above.
  text_size_mm <- function(multiplier = 1) {
    base_font_size * multiplier / 2.845276
  }

  title_size_mm <- function(multiplier = 1) {
    title_font_size * multiplier / 2.845276
  }

  panel_title_size_mm <- function(multiplier = 1) {
    panel_title_font_size * multiplier / 2.845276
  }

  wrap_label <- function(values, width = 18L) {
    vapply(
      values,
      function(value) {
        paste(strwrap(as.character(value), width = width), collapse = "\n")
      },
      character(1)
    )
  }

  format_trait_label_a <- function(values) {
    labels <- gsub("_", " ", as.character(values), fixed = TRUE)
    labels[labels == "Leaf Vcmax per dry mass"] <-
      "Leaf Vcmax per\ndry mass"
    labels
  }

  format_trait_label_b <- function(values) {
    labels <- gsub("_", " ", as.character(values), fixed = TRUE)
    labels[labels == "Stem conduit diameter"] <-
      "Stem conduit\ndiameter"
    labels[labels == "Leaf Vcmax per dry mass"] <-
      "Leaf Vcmax per\ndry mass"
    labels[labels == "Stomatal conductance"] <-
      "Stomatal\nconductance"
    labels
  }

  if (
    is.character(r2_data) &&
    length(r2_data) == 1L &&
    !is.na(r2_data)
  ) {
    if (!file.exists(r2_data)) {
      stop("The R2 data file does not exist.", call. = FALSE)
    }

    r2_data <- utils::read.csv(
      r2_data,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  if (
    !is.numeric(title_font_size) ||
    length(title_font_size) != 1L ||
    is.na(title_font_size) ||
    title_font_size <= 0
  ) {
    stop(
      "title_font_size must be one positive number.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(panel_title_font_size) ||
    length(panel_title_font_size) != 1L ||
    is.na(panel_title_font_size) ||
    panel_title_font_size <= 0
  ) {
    stop(
      "panel_title_font_size must be one positive number.",
      call. = FALSE
    )
  }

  if (!is.data.frame(r2_data)) {
    stop(
      "r2_data must be a data frame or a CSV file path.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(trait_spacing_multiplier) ||
    length(trait_spacing_multiplier) != 1L ||
    is.na(trait_spacing_multiplier) ||
    trait_spacing_multiplier <= 0
  ) {
    stop(
      "trait_spacing_multiplier must be one positive number.",
      call. = FALSE
    )
  }

  required_columns <- c(
    "group",
    "trait",
    "analysis",
    "configuration_label",
    "folds_used",
    "r2"
  )

  missing_columns <- setdiff(required_columns, names(r2_data))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  group_order <- c(
    "All8",
    paste("Cluster", 1:6)
  )

  group_display <- c(
    All8 = "(I) 8 selected traits",
    `Cluster 1` = "(II) Cluster I",
    `Cluster 2` = "(III) Cluster II",
    `Cluster 3` = "(IV) Cluster III",
    `Cluster 4` = "(V) Cluster IV",
    `Cluster 5` = "(VI) Cluster V",
    `Cluster 6` = "(VII) Cluster VI"
  )

  group_colors <- c(
    All8 = "#4D4D4D",
    `Cluster 1` = "#D89000",
    `Cluster 2` = "#E64291",
    `Cluster 3` = "#7389C6",
    `Cluster 4` = "#4C9B35",
    `Cluster 5` = "#3579B9",
    `Cluster 6` = "#149977"
  )

  trait_order <- list(
    All8 = c(
      "Wood_density",
      "Specific_leaf_area",
      "Seed_dry_mass",
      "Leaf_P_per_mass",
      "Stem_conduit_diameter",
      "Tree_height",
      "Root_depth",
      "Bark_thickness"
    ),
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

  configurations <- c(
    "BM_FTF",
    "BM_TTF",
    "lambda_FTF",
    "lambda_TTF"
  )

  configuration_display <- c(
    BM_FTF = "BM · FTF",
    BM_TTF = "BM · TTF",
    lambda_FTF = "lambda · FTF",
    lambda_TTF = "lambda · TTF"
  )

  expected_traits <- unname(unlist(trait_order, use.names = FALSE))
  expected_grid <- do.call(
    rbind,
    lapply(
      group_order,
      function(group_name) {
        expand.grid(
          group = group_name,
          trait = trait_order[[group_name]],
          configuration_label = configurations,
          analysis = c("Single trait", "Multi-trait"),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  observed_key <- paste(
    r2_data$group,
    r2_data$trait,
    r2_data$configuration_label,
    r2_data$analysis,
    sep = "\r"
  )
  expected_key <- paste(
    expected_grid$group,
    expected_grid$trait,
    expected_grid$configuration_label,
    expected_grid$analysis,
    sep = "\r"
  )

  valid_analysis <- r2_data$analysis %in% c(
    "Single trait",
    "Multi-trait"
  )
  valid_values <- is.numeric(r2_data$r2) &&
    !anyNA(r2_data$r2) &&
    all(is.finite(r2_data$r2))

  if (
    !all(valid_analysis) ||
    !valid_values ||
    anyDuplicated(observed_key) ||
    !setequal(observed_key, expected_key)
  ) {
    stop(
      paste0(
        "The R2 data must contain exactly one finite single-trait and ",
        "one finite multi-trait result for every group, trait, and ",
        "configuration."
      ),
      call. = FALSE
    )
  }

  single_data <- r2_data[
    r2_data$analysis == "Single trait",
    c(
      "group",
      "trait",
      "configuration_label",
      "folds_used",
      "r2"
    ),
    drop = FALSE
  ]
  names(single_data)[names(single_data) == "r2"] <- "single_r2"

  multi_data <- r2_data[
    r2_data$analysis == "Multi-trait",
    c(
      "group",
      "trait",
      "configuration_label",
      "folds_used",
      "r2"
    ),
    drop = FALSE
  ]
  names(multi_data)[names(multi_data) == "r2"] <- "multi_r2"

  comparison_data <- merge(
    single_data,
    multi_data,
    by = c(
      "group",
      "trait",
      "configuration_label",
      "folds_used"
    ),
    all = TRUE,
    sort = FALSE
  )

  if (
    nrow(comparison_data) != nrow(expected_grid) / 2L ||
    anyNA(comparison_data$single_r2) ||
    anyNA(comparison_data$multi_r2)
  ) {
    stop(
      "Single- and multi-trait results could not be paired.",
      call. = FALSE
    )
  }

  comparison_data$delta_r2 <- (
    comparison_data$multi_r2 - comparison_data$single_r2
  )
  comparison_data$direction <- ifelse(
    comparison_data$delta_r2 > 1e-12,
    "Multi higher",
    ifelse(
      comparison_data$delta_r2 < -1e-12,
      "Single higher",
      "Equal"
    )
  )
  comparison_data$model <- ifelse(
    grepl("^BM_", comparison_data$configuration_label),
    "BM",
    "lambda"
  )
  comparison_data$pheno_error <- ifelse(
    grepl("_TTF$", comparison_data$configuration_label),
    "TRUE",
    "FALSE"
  )
  comparison_data$group_display <- unname(
    group_display[comparison_data$group]
  )
  comparison_data$trait_display <- gsub(
    "_",
    " ",
    comparison_data$trait,
    fixed = TRUE
  )
  comparison_data$configuration_display <- unname(
    configuration_display[comparison_data$configuration_label]
  )

  group_rank <- match(comparison_data$group, group_order)
  trait_rank <- vapply(
    seq_len(nrow(comparison_data)),
    function(index) {
      match(
        comparison_data$trait[[index]],
        trait_order[[comparison_data$group[[index]]]]
      )
    },
    integer(1)
  )
  configuration_rank <- match(
    comparison_data$configuration_label,
    configurations
  )

  comparison_data <- comparison_data[
    order(group_rank, trait_rank, configuration_rank),
    ,
    drop = FALSE
  ]
  rownames(comparison_data) <- NULL

  compress_negative <- function(values) {
    ifelse(
      values < 0,
      values * negative_compression,
      values
    )
  }

  comparison_data$single_plot <- compress_negative(
    comparison_data$single_r2
  )
  comparison_data$multi_plot <- compress_negative(
    comparison_data$multi_r2
  )

  positive_maximum <- max(
    0.1,
    comparison_data$single_r2[comparison_data$single_r2 >= 0],
    comparison_data$multi_r2[comparison_data$multi_r2 >= 0],
    na.rm = TRUE
  )
  axis_maximum <- ceiling(positive_maximum * 10) / 10
  positive_breaks <- seq(0, axis_maximum, by = 0.1)
  negative_minimum <- min(
    0,
    comparison_data$single_plot,
    comparison_data$multi_plot,
    na.rm = TRUE
  )
  axis_minimum <- negative_minimum * 1.30
  r2_label_x <- axis_minimum * 0.75
  axis_breaks <- c(r2_label_x, positive_breaks)
  axis_labels <- c(
    "R²",
    formatC(positive_breaks, format = "f", digits = 1)
  )
  delta_x <- axis_maximum + 0.010
  plotting_maximum <- axis_maximum + 0.17

  comparison_data$delta_label <- sprintf(
    "%+.3f",
    comparison_data$delta_r2
  )
  comparison_data$single_label <- ifelse(
    comparison_data$single_r2 < 0,
    sprintf("%.2f", comparison_data$single_r2),
    ""
  )
  comparison_data$multi_label <- ifelse(
    comparison_data$multi_r2 < 0,
    sprintf("%.2f", comparison_data$multi_r2),
    ""
  )

  direction_colors <- c(
    `Multi higher` = "#2878B5",
    `Single higher` = "#D55E00",
    Equal = "#6B6B6B"
  )

  configuration_shapes <- c(
    BM_FTF = 21,
    BM_TTF = 24,
    lambda_FTF = 22,
    lambda_TTF = 23
  )

  unique_traits <- unique(unname(unlist(trait_order, use.names = FALSE)))
  trait_colors <- stats::setNames(
    grDevices::hcl.colors(length(unique_traits), palette = "Dark 3"),
    unique_traits
  )

  make_dumbbell_group <- function(group_name) {
    current <- comparison_data[
      comparison_data$group == group_name,
      ,
      drop = FALSE
    ]
    current$trait_index <- match(
      current$trait,
      trait_order[[group_name]]
    )
    current$configuration_index <- match(
      current$configuration_label,
      configurations
    )
    current$negative_value <- pmin(
      current$single_r2,
      current$multi_r2
    )
    current$negative_label <- ifelse(
      current$negative_value < 0,
      sprintf("%.2f", current$negative_value),
      ""
    )
    current$negative_plot <- pmin(
      current$single_plot,
      current$multi_plot
    )
    number_of_traits <- length(trait_order[[group_name]])
    block_height <- 5
    current$y <- (
      (number_of_traits - current$trait_index) * block_height +
        (5 - current$configuration_index)
    )

    trait_table <- data.frame(
      trait = trait_order[[group_name]],
      trait_index = seq_len(number_of_traits),
      stringsAsFactors = FALSE
    )
    trait_table$label <- format_trait_label_a(trait_table$trait)
    trait_table$ymin <- (
      (number_of_traits - trait_table$trait_index) * block_height + 0.45
    )
    trait_table$ymax <- (
      (number_of_traits - trait_table$trait_index) * block_height + 4.55
    )
    trait_table$ymid <- (trait_table$ymin + trait_table$ymax) / 2
    trait_table$label_color <- unname(trait_colors[trait_table$trait])
    separator_data <- data.frame(
      y = if (number_of_traits > 1L) {
        rev(seq_len(number_of_traits - 1L)) * block_height
      } else {
        numeric(0)
      }
    )

    data_span <- axis_maximum - axis_minimum
    current$negative_label_x <- (
      current$negative_plot - data_span * 0.020
    )
    table_gap <- data_span * 0.035
    table_width <- data_span * 0.64
    table_xmax <- axis_minimum - table_gap
    table_xmin <- table_xmax - table_width
    panel_xmax <- axis_maximum + data_span * 0.035
    panel_ymax <- number_of_traits * block_height + 0.55
    show_shared_axis <- group_name %in% c("All8", "Cluster 1")

    ggplot2::ggplot() +
      ggplot2::geom_rect(
        data = trait_table,
        ggplot2::aes(
          xmin = table_xmin,
          xmax = table_xmax,
          ymin = ymin,
          ymax = ymax
        ),
        inherit.aes = FALSE,
        fill = "white",
        color = NA,
        linewidth = 0
      ) +
      ggplot2::geom_segment(
        data = separator_data,
        ggplot2::aes(
          x = table_xmin,
          xend = panel_xmax,
          y = y,
          yend = y
        ),
        inherit.aes = FALSE,
        color = "grey82",
        linewidth = 0.28
      ) +
      ggplot2::geom_text(
        data = trait_table,
        ggplot2::aes(
          x = table_xmin + table_width * 0.02,
          y = ymid,
          label = label,
          color = label_color
        ),
        inherit.aes = FALSE,
        hjust = 0,
        lineheight = 0.92,
        size = title_size_mm(),
        family = font_family,
        show.legend = FALSE
      ) +
      ggplot2::geom_segment(
        data = data.frame(x = positive_breaks),
        ggplot2::aes(
          x = x,
          xend = x,
          y = 0.45,
          yend = panel_ymax
        ),
        inherit.aes = FALSE,
        color = "grey92",
        linewidth = 0.18
      ) +
      ggplot2::annotate(
        "segment",
        x = 0,
        xend = 0,
        y = 0.45,
        yend = panel_ymax,
        color = "grey48",
        linewidth = 0.30
      ) +
      ggplot2::geom_segment(
        data = current,
        ggplot2::aes(
          x = single_plot,
          xend = multi_plot,
          y = y,
          yend = y,
          color = direction
        ),
        inherit.aes = FALSE,
        linewidth = 0.30,
        alpha = 0.95
      ) +
      ggplot2::geom_point(
        data = current,
        ggplot2::aes(
          x = single_plot,
          y = y,
          shape = configuration_label,
          color = direction
        ),
        inherit.aes = FALSE,
        size = 1.52,
        stroke = 0.42,
        fill = "white",
        show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = current,
        ggplot2::aes(
          x = multi_plot,
          y = y,
          shape = configuration_label,
          fill = direction,
          color = direction
        ),
        inherit.aes = FALSE,
        size = 1.62,
        stroke = 0.42,
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = current[current$negative_value < 0, , drop = FALSE],
        ggplot2::aes(
          x = negative_label_x,
          y = y,
          label = negative_label,
          color = direction
        ),
        inherit.aes = FALSE,
        hjust = 1,
        vjust = 0.5,
        size = text_size_mm(0.82),
        family = font_family,
        show.legend = FALSE
      ) +
      ggplot2::scale_color_manual(
        values = c(
          direction_colors,
          stats::setNames(trait_colors, trait_colors)
        ),
        guide = "none"
      ) +
      ggplot2::scale_fill_manual(
        values = direction_colors,
        guide = "none"
      ) +
      ggplot2::scale_shape_manual(
        values = configuration_shapes,
        guide = "none"
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0.35, panel_ymax),
        expand = c(0, 0)
      ) +
      ggplot2::scale_x_continuous(
        breaks = axis_breaks,
        labels = axis_labels,
        limits = c(table_xmin, panel_xmax),
        position = "top",
        expand = c(0, 0)
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(
        title = unname(group_display[[group_name]]),
        x = NULL,
        y = NULL
      ) +
      ggplot2::theme_minimal(
        base_size = base_font_size,
        base_family = font_family
      ) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.x = if (show_shared_axis) {
          ggplot2::element_text(
            color = "grey40",
            size = base_font_size
          )
        } else {
          ggplot2::element_blank()
        },
        axis.ticks.x = if (show_shared_axis) {
          ggplot2::element_line(
            color = "grey55",
            linewidth = 0.25
          )
        } else {
          ggplot2::element_blank()
        },
        axis.title.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          face = "plain",
          family = font_family,
          color = unname(group_colors[[group_name]]),
          size = title_font_size,
          hjust = 0
        ),
        plot.background = ggplot2::element_rect(
          fill = NA,
          color = unname(group_colors[[group_name]]),
          linewidth = 0.75
        ),
        legend.position = "none",
        plot.margin = ggplot2::margin(3, 4, 3, 3, unit = "pt")
      )
  }

  group_plots <- lapply(group_order, make_dumbbell_group)
  names(group_plots) <- group_order

  shared_axis_groups <- c("All8", "Cluster 1")
  group_height_units <- vapply(
    group_order,
    function(group_name) {
      length(trait_order[[group_name]]) * trait_spacing_multiplier +
        0.70 +
        ifelse(group_name %in% shared_axis_groups, 0.20, 0)
    },
    numeric(1)
  )
  cluster_gap_unit <- 0.10
  left_groups <- c("All8", "Cluster 2")
  right_groups <- c("Cluster 1", paste("Cluster", 3:6))
  left_total <- sum(group_height_units[left_groups]) + cluster_gap_unit
  right_total <- sum(group_height_units[right_groups]) +
    (length(right_groups) - 1L) * cluster_gap_unit
  column_total <- max(left_total, right_total)

  left_column <- patchwork::wrap_plots(
    group_plots[["All8"]],
    patchwork::plot_spacer(),
    group_plots[["Cluster 2"]],
    patchwork::plot_spacer(),
    ncol = 1,
    heights = c(
      group_height_units[["All8"]],
      cluster_gap_unit,
      group_height_units[["Cluster 2"]],
      max(0.001, column_total - left_total)
    )
  )
  right_column <- patchwork::wrap_plots(
    group_plots[["Cluster 1"]],
    patchwork::plot_spacer(),
    group_plots[["Cluster 3"]],
    patchwork::plot_spacer(),
    group_plots[["Cluster 4"]],
    patchwork::plot_spacer(),
    group_plots[["Cluster 5"]],
    patchwork::plot_spacer(),
    group_plots[["Cluster 6"]],
    patchwork::plot_spacer(),
    ncol = 1,
    heights = c(
      group_height_units[["Cluster 1"]],
      cluster_gap_unit,
      group_height_units[["Cluster 3"]],
      cluster_gap_unit,
      group_height_units[["Cluster 4"]],
      cluster_gap_unit,
      group_height_units[["Cluster 5"]],
      cluster_gap_unit,
      group_height_units[["Cluster 6"]],
      max(0.001, column_total - right_total)
    )
  )

  configuration_key <- data.frame(
    x = c(0.20, 0.36, 0.20, 0.36),
    label_x = c(0.225, 0.385, 0.225, 0.385),
    y = c(0.66, 0.66, 0.34, 0.34),
    configuration = configurations,
    label = unname(configuration_display[configurations]),
    stringsAsFactors = FALSE
  )

  legend_panel <- ggplot2::ggplot() +
    ggplot2::annotate(
      "rect", xmin = 0.025, xmax = 0.975, ymin = 0.08, ymax = 0.92,
      fill = "white", color = "grey72", linewidth = 0.45
    ) +
    ggplot2::annotate(
      "segment", x = c(0.16, 0.50, 0.68), xend = c(0.16, 0.50, 0.68),
      y = 0.12, yend = 0.88, color = "grey84", linewidth = 0.30
    ) +
    ggplot2::annotate(
      "text", x = 0.0925, y = 0.52, label = "Context shape",
      hjust = 0.5, size = text_size_mm(), color = "grey25"
    ) +
    ggplot2::geom_point(
      data = configuration_key,
      ggplot2::aes(x = x, y = y, shape = configuration),
      inherit.aes = FALSE, size = 2.15, stroke = 0.50,
      fill = "white", color = "grey20"
    ) +
    ggplot2::geom_text(
      data = configuration_key,
      ggplot2::aes(x = label_x, y = y, label = label),
      inherit.aes = FALSE, hjust = 0, size = text_size_mm(),
      family = font_family, color = "grey25"
    ) +
    ggplot2::annotate(
      "text", x = 0.59, y = 0.70, label = "Endpoint fill",
      hjust = 0.5, size = text_size_mm(), color = "grey25"
    ) +
    ggplot2::annotate(
      "text", x = 0.59, y = 0.28,
      label = expression(Delta * R^2 * " direction"),
      hjust = 0.5, size = text_size_mm(), color = "grey25"
    ) +
    ggplot2::annotate(
      "point", x = 0.715, y = 0.70, shape = 21, size = 2.15,
      stroke = 0.50, fill = "white", color = "grey20"
    ) +
    ggplot2::annotate(
      "text", x = 0.738, y = 0.70, label = "Univariate",
      hjust = 0, size = text_size_mm(0.90), color = "grey25"
    ) +
    ggplot2::annotate(
      "point", x = 0.855, y = 0.70, shape = 21, size = 2.15,
      stroke = 0.50, fill = "grey20", color = "grey20"
    ) +
    ggplot2::annotate(
      "text", x = 0.878, y = 0.70, label = "Multivariate",
      hjust = 0, size = text_size_mm(0.90), color = "grey25"
    ) +
    ggplot2::annotate(
      "segment", x = 0.705, xend = 0.735, y = 0.28, yend = 0.28,
      linewidth = 0.42, color = direction_colors[["Multi higher"]]
    ) +
    ggplot2::annotate(
      "text", x = 0.742, y = 0.28, label = "Multivariate\nbetter",
      hjust = 0, lineheight = 0.85,
      size = text_size_mm(0.90), color = "grey25"
    ) +
    ggplot2::annotate(
      "segment", x = 0.842, xend = 0.872, y = 0.28, yend = 0.28,
      linewidth = 0.42, color = direction_colors[["Single higher"]]
    ) +
    ggplot2::annotate(
      "text", x = 0.879, y = 0.28, label = "Univariate\nbetter",
      hjust = 0, lineheight = 0.85,
      size = text_size_mm(0.90), color = "grey25"
    ) +
    ggplot2::scale_shape_manual(values = configuration_shapes, guide = "none") +
    ggplot2::coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1),
      clip = "off"
    ) +
    ggplot2::theme_void()

  dumbbell_plot <- (
    patchwork::wrap_plots(
      left_column,
      patchwork::plot_spacer(),
      right_column,
      ncol = 3,
      widths = c(1.00, 0.035, 1.00)
    ) /
      legend_panel
  ) +
    patchwork::plot_layout(heights = c(0.87, 0.13)) +
    patchwork::plot_annotation(
      title = "Univariate versus multivariate performance (R² and ΔR²)",
      subtitle = NULL,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "plain",
          size = title_font_size,
          hjust = 0
        ),
        plot.subtitle = ggplot2::element_text(
          color = "grey38",
          size = base_font_size,
          hjust = 0,
          margin = ggplot2::margin(b = 8)
        )
      )
    )

  if (is.null(delta_fill_limit)) {
    delta_fill_limit <- max(
      0.05,
      stats::quantile(
        abs(comparison_data$delta_r2),
        probs = 0.90,
        na.rm = TRUE,
        names = FALSE
      )
    )
  }

  comparison_data$delta_fill <- pmax(
    -delta_fill_limit,
    pmin(delta_fill_limit, comparison_data$delta_r2)
  )
  comparison_data$heatmap_delta_label <- sprintf(
    "%+.3f",
    comparison_data$delta_r2
  )
  comparison_data$model <- factor(
    comparison_data$model,
    levels = c("lambda", "BM")
  )
  comparison_data$pheno_error <- factor(
    comparison_data$pheno_error,
    levels = c("FALSE", "TRUE"),
    labels = c("FTF", "TTF")
  )

  make_heatmap_group <- function(group_name) {
    current <- comparison_data[
      comparison_data$group == group_name,
      ,
      drop = FALSE
    ]
    current$trait_display <- factor(
      current$trait_display,
      levels = gsub(
        "_",
        " ",
        trait_order[[group_name]],
        fixed = TRUE
      )
    )
    levels(current$trait_display) <- format_trait_label_b(
      levels(current$trait_display)
    )
    trait_border_data <- unique(
      current[c("trait", "trait_display")]
    )
    trait_border_data$trait_color <- unname(
      trait_colors[trait_border_data$trait]
    )
    trait_border_data$cell_xmin <- 0.50
    trait_border_data$cell_xmax <- 2.50
    trait_border_data$cell_ymin <- 0.50
    trait_border_data$cell_ymax <- 2.50
    trait_border_data$header_ymin <- 2.55
    trait_border_data$header_ymax <- 4.05
    trait_border_data$header_x <- 1.50
    trait_border_data$header_y <- 3.30
    current$heatmap_x <- as.numeric(current$pheno_error)
    current$heatmap_y <- as.numeric(current$model)

    ggplot2::ggplot(
      current,
      ggplot2::aes(
        x = heatmap_x,
        y = heatmap_y,
        fill = delta_fill
      )
    ) +
      ggplot2::geom_tile(
        color = "white",
        linewidth = 0.18,
        width = 0.985,
        height = 0.92
      ) +
      ggplot2::geom_rect(
        data = trait_border_data,
        ggplot2::aes(
          xmin = cell_xmin,
          xmax = cell_xmax,
          ymin = cell_ymin,
          ymax = cell_ymax,
          color = trait_color
        ),
        inherit.aes = FALSE,
        fill = NA,
        linewidth = 0.55
      ) +
      ggplot2::geom_rect(
        data = trait_border_data,
        ggplot2::aes(
          xmin = cell_xmin,
          xmax = cell_xmax,
          ymin = header_ymin,
          ymax = header_ymax
        ),
        inherit.aes = FALSE,
        fill = "grey96",
        color = NA
      ) +
      ggplot2::geom_text(
        data = trait_border_data,
        ggplot2::aes(
          x = header_x,
          y = header_y,
          label = trait_display,
          color = trait_color
        ),
        inherit.aes = FALSE,
        size = title_size_mm(0.84),
        lineheight = 0.90,
        family = font_family,
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = heatmap_delta_label),
        size = text_size_mm(),
        fontface = "plain",
        color = "grey5"
      ) +
      ggplot2::facet_wrap(
        ggplot2::vars(trait_display),
        ncol = min(4L, length(trait_order[[group_name]])),
        scales = "fixed"
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0.34, 2.66),
        expand = c(0, 0)
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0.30, 4.25),
        expand = c(0, 0)
      ) +
      ggplot2::scale_fill_gradient2(
        low = "#C84A00",
        mid = "#F7F7F7",
        high = "#1769AA",
        midpoint = 0,
        limits = c(-delta_fill_limit, delta_fill_limit),
        breaks = c(
          -delta_fill_limit,
          -delta_fill_limit / 2,
          0,
          delta_fill_limit / 2,
          delta_fill_limit
        ),
        oob = scales::squish,
        name = expression(Delta * R^2)
      ) +
      ggplot2::scale_color_identity(guide = "none") +
      ggplot2::labs(
        title = unname(group_display[[group_name]]),
        x = NULL,
        y = NULL
      ) +
      ggplot2::theme_minimal(
        base_size = base_font_size,
        base_family = font_family
      ) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        strip.text = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        panel.spacing.x = grid::unit(1.2, "pt"),
        panel.spacing.y = grid::unit(2.0, "pt"),
        plot.title = ggplot2::element_text(
          face = "plain",
          family = font_family,
          color = unname(group_colors[[group_name]]),
          size = title_font_size,
          hjust = 0,
          margin = ggplot2::margin(
            t = 4,
            r = 0,
            b = 3,
            l = 4,
            unit = "pt"
          )
        ),
        plot.background = ggplot2::element_rect(
          fill = NA,
          color = unname(group_colors[[group_name]]),
          linewidth = 0.75
        ),
        legend.position = "bottom",
        legend.key.width = grid::unit(2.0, "cm"),
        plot.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt")
      )
  }

  group_heatmaps <- lapply(group_order, make_heatmap_group)
  names(group_heatmaps) <- group_order

  heatmap_group_stack <- patchwork::wrap_plots(
    group_heatmaps,
    ncol = 1,
    heights = vapply(
      group_order,
      function(group_name) {
        max(1.0, length(trait_order[[group_name]]) / 4)
      },
      numeric(1)
    ),
    guides = "collect"
  )

  heatmap_key_data <- expand.grid(
    pheno_error = factor(
      c("FTF", "TTF"),
      levels = c("FTF", "TTF")
    ),
    model = factor(
      c("lambda", "BM"),
      levels = c("lambda", "BM")
    )
  )

  heatmap_layout_key <- ggplot2::ggplot(
    heatmap_key_data,
    ggplot2::aes(x = pheno_error, y = model)
  ) +
    ggplot2::geom_tile(
      fill = "white",
      color = "grey78",
      linewidth = 0.24,
      width = 0.97,
      height = 0.94
    ) +
    ggplot2::labs(
      title = "Heatmap layout",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = base_font_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(
        color = "grey38",
        size = base_font_size
      ),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "plain",
        color = "grey28",
        size = base_font_size,
        hjust = 0.5,
        margin = ggplot2::margin(b = 2, unit = "pt")
      ),
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA,
        linewidth = 0.45
      ),
      plot.margin = ggplot2::margin(1, 1, 1, 1, unit = "pt")
    )

  heatmap_key_row <- patchwork::wrap_plots(
    patchwork::plot_spacer(),
    heatmap_layout_key,
    patchwork::plot_spacer(),
    ncol = 3,
    widths = c(1, 0.42, 1)
  )

  heatmap_color_key_data <- data.frame(
    delta = seq(
      -delta_fill_limit,
      delta_fill_limit,
      length.out = 201L
    ),
    y = 1
  )

  heatmap_color_key <- ggplot2::ggplot(
    heatmap_color_key_data,
    ggplot2::aes(x = delta, y = y, fill = delta)
  ) +
    ggplot2::geom_tile(
      width = 2 * delta_fill_limit / 200,
      height = 0.12
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#C84A00",
      mid = "#F7F7F7",
      high = "#1769AA",
      midpoint = 0,
      limits = c(-delta_fill_limit, delta_fill_limit),
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(
        -delta_fill_limit,
        -delta_fill_limit / 2,
        0,
        delta_fill_limit / 2,
        delta_fill_limit
      ),
      labels = function(values) {
        formatC(values, format = "f", digits = 3)
      },
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(ylim = c(0.90, 1.10)) +
    ggplot2::labs(
      title = "Delta R² scale",
      x = expression(Delta * R^2),
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = base_font_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        color = "grey38",
        size = base_font_size,
        margin = ggplot2::margin(t = 1, unit = "pt")
      ),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(
        color = "grey28",
        size = base_font_size,
        margin = ggplot2::margin(t = 1, unit = "pt")
      ),
      plot.title = ggplot2::element_text(
        face = "plain",
        color = "grey28",
        size = base_font_size,
        hjust = 0.5,
        margin = ggplot2::margin(b = 2, unit = "pt")
      ),
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA,
        linewidth = 0.45
      ),
      plot.margin = ggplot2::margin(5, 5, 5, 5, unit = "pt")
    )

  heatmap_guide_row <- patchwork::wrap_plots(
    patchwork::plot_spacer(),
    heatmap_color_key,
    patchwork::plot_spacer(),
    heatmap_layout_key,
    patchwork::plot_spacer(),
    ncol = 5,
    widths = c(0.18, 0.72, 0.18, 0.58, 0.34)
  )

  heatmap_plot <- patchwork::wrap_plots(
    heatmap_group_stack,
    heatmap_key_row,
    ncol = 1,
    heights = c(10, 1.15),
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = expression(Delta * R^2 * " summary: multivariate minus univariate"),
      subtitle = paste0(
        "Cells show ΔR² only; color is clipped symmetrically at ±",
        formatC(delta_fill_limit, format = "f", digits = 3),
        "."
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "plain",
          size = title_font_size,
          hjust = 0
        ),
        plot.subtitle = ggplot2::element_text(
          color = "grey38",
          size = base_font_size,
          hjust = 0,
          margin = ggplot2::margin(b = 8)
        )
      )
    ) &
    ggplot2::theme(legend.position = "bottom")

  make_panel_label <- function(label) {
    ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0.5,
        label = label,
        hjust = 0,
        vjust = 0.5,
        fontface = "plain",
        family = font_family,
        size = panel_title_size_mm(),
        color = "grey8"
      ) +
      ggplot2::coord_cartesian(
        xlim = c(0, 1),
        ylim = c(0, 1),
        clip = "off"
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0, 3, 0, 3, unit = "pt")
      )
  }

  dumbbell_body <- patchwork::wrap_plots(
    left_column,
    patchwork::plot_spacer(),
    right_column,
    ncol = 3,
    widths = c(1.00, 0.035, 1.00)
  )

  combined_panel_a <- (
    make_panel_label(
      "A  Univariate versus multivariate R² performance"
    ) /
      dumbbell_body /
      legend_panel
  ) +
    patchwork::plot_layout(
      heights = c(
        1.35,
        18.2 * trait_spacing_multiplier,
        2.30
      )
    )

  all8_heatmap_combined <- group_heatmaps[["All8"]] +
    ggplot2::theme(legend.position = "none")

  cluster_heatmap_plots <- lapply(
    group_heatmaps[paste("Cluster", 1:6)],
    function(plot) {
      plot + ggplot2::theme(legend.position = "none")
    }
  )

  cluster_heatmap_grid <- patchwork::wrap_plots(
    cluster_heatmap_plots,
    ncol = 2,
    heights = c(1, 1, 1),
    guides = "keep"
  )

  combined_panel_b <- patchwork::wrap_plots(
    make_panel_label(
      "B  ΔR² summary (ΔR² = multivariate R² − univariate R²)"
    ),
    all8_heatmap_combined,
    patchwork::plot_spacer(),
    cluster_heatmap_grid,
    heatmap_guide_row,
    ncol = 1,
    heights = c(0.68, 4.1, 0.12, 8.2, 1.15),
    guides = "keep"
  )

  combined_plot <- patchwork::wrap_plots(
    combined_panel_a,
    patchwork::plot_spacer(),
    combined_panel_b,
    ncol = 1,
    heights = c(6.25, 0.12, 6.63)
  )

  save_plot <- function(plot, path, save_width, save_height) {
    if (is.null(path)) {
      return(invisible(NULL))
    }

    if (
      !is.character(path) ||
      length(path) != 1L ||
      is.na(path) ||
      !nzchar(trimws(path))
    ) {
      stop("Output paths must be non-empty strings.", call. = FALSE)
    }

    extension <- tolower(tools::file_ext(path))
    if (!extension %in% c("png", "pdf")) {
      stop("Plot output must be PNG or PDF.", call. = FALSE)
    }

    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    arguments <- list(
      filename = path,
      plot = plot,
      width = save_width,
      height = save_height,
      units = "in",
      dpi = dpi,
      bg = "white"
    )

    if (
      extension == "png" &&
      requireNamespace("ragg", quietly = TRUE)
    ) {
      arguments$device <- ragg::agg_png
    }

    do.call(ggplot2::ggsave, arguments)
    invisible(path)
  }

  detailed_table <- comparison_data[
    ,
    c(
      "group_display",
      "trait",
      "model",
      "pheno_error",
      "configuration_label",
      "folds_used",
      "single_r2",
      "multi_r2",
      "delta_r2",
      "direction"
    ),
    drop = FALSE
  ]
  names(detailed_table)[1L] <- "group"

  summary_rows <- lapply(
    split(
      seq_len(nrow(comparison_data)),
      interaction(
        comparison_data$group,
        comparison_data$trait,
        drop = TRUE,
        lex.order = TRUE
      )
    ),
    function(indices) {
      current <- comparison_data[indices, , drop = FALSE]
      current <- current[
        match(configurations, current$configuration_label),
        ,
        drop = FALSE
      ]

      data.frame(
        group = current$group_display[[1L]],
        trait = current$trait[[1L]],
        BM_FTF_delta = current$delta_r2[[1L]],
        BM_TTF_delta = current$delta_r2[[2L]],
        lambda_FTF_delta = current$delta_r2[[3L]],
        lambda_TTF_delta = current$delta_r2[[4L]],
        mean_delta = mean(current$delta_r2),
        multi_higher_count = sum(current$delta_r2 > 1e-12),
        single_higher_count = sum(current$delta_r2 < -1e-12),
        stringsAsFactors = FALSE
      )
    }
  )
  summary_table <- do.call(rbind, summary_rows)
  rownames(summary_table) <- NULL

  summary_group_rank <- match(
    summary_table$group,
    unname(group_display[group_order])
  )
  summary_trait_rank <- vapply(
    seq_len(nrow(summary_table)),
    function(index) {
      original_group <- group_order[[summary_group_rank[[index]]]]
      match(summary_table$trait[[index]], trait_order[[original_group]])
    },
    integer(1)
  )
  summary_table <- summary_table[
    order(summary_group_rank, summary_trait_rank),
    ,
    drop = FALSE
  ]

  write_table <- function(table, path) {
    if (is.null(path)) {
      return(invisible(NULL))
    }

    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(table, path, row.names = FALSE, na = "")
    invisible(path)
  }

  save_plot(dumbbell_plot, output_file, width, height)
  save_plot(
    heatmap_plot,
    heatmap_file,
    width,
    18
  )
  save_plot(
    combined_plot,
    combined_file,
    combined_width,
    combined_height
  )
  write_table(detailed_table, detailed_table_file)
  write_table(summary_table, summary_table_file)

  if (interactive()) {
    print(combined_plot)
  }

  invisible(
    list(
      dumbbell_plot = dumbbell_plot,
      heatmap_plot = heatmap_plot,
      combined_plot = combined_plot,
      detailed_table = detailed_table,
      summary_table = summary_table,
      comparison_data = comparison_data
    )
  )
}

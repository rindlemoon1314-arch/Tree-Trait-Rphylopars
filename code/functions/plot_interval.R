plot_single_lambda_ttf_interval_calibration <- function(
    metrics_input = NULL,
    single_dir = NULL,
    calculate_metrics_file = paste0(
      here::here("code/functions/"),
      "calculate_phy_metrics.R"
    ),
    folds = c(1L, 3L, 5L, 6L, 8L, 10L),
    traits = NULL,
    interval = 95,
    use_transformed_truth = TRUE,
    traits_per_page = 18L,
    max_points_per_trait = Inf,
    display_seed = 123L,
    output_dir = NULL,
    file_prefix = "single_lambda_TTF_interval_calibration",
    width = 8.25,
    height = 10.8,
    dpi = 300,
    base_size = 9.5,
    horizontal_padding = 0.005
) {
  # Use the same Word-compatible page width as the other manuscript figures
  # and keep the total figure height below the 11.25-inch page limit.
  width <- 8.25
  height <- min(as.numeric(height)[1L], 11.25)
  if (!is.finite(height) || height <= 0) {
    stop("height must be a positive finite number.", call. = FALSE)
  }
  horizontal_padding <- as.numeric(horizontal_padding)[1L]
  if (!is.finite(horizontal_padding) || horizontal_padding < 0 ||
      horizontal_padding >= 0.25) {
    stop(
      "horizontal_padding must be between 0 and 0.25.",
      call. = FALSE
    )
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.", call. = FALSE)
  }

  if (!is.numeric(interval) || length(interval) != 1L ||
      is.na(interval) || interval <= 0 || interval >= 100) {
    stop("interval must be one percentage strictly between 0 and 100.", call. = FALSE)
  }

  traits_per_page <- as.integer(traits_per_page)
  if (is.na(traits_per_page) || traits_per_page < 1L) {
    stop("traits_per_page must be a positive integer.", call. = FALSE)
  }

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Use an existing metrics object when supplied. Otherwise calculate the
  # species-level assessment directly from all saved single-trait lambda-TTF
  # PHY objects. calculate_phy_metrics() averages held-out records within each
  # fold, species, and trait before calculating intervals and coverage.
  if (!is.null(metrics_input)) {
    metrics <- if (is.character(metrics_input) && length(metrics_input) == 1L) {
      if (!file.exists(metrics_input)) {
        stop("metrics_input does not exist.", call. = FALSE)
      }
      readRDS(metrics_input)
    } else {
      metrics_input
    }
  } else {
    if (is.null(single_dir) || length(single_dir) < 1L ||
        any(!dir.exists(single_dir))) {
      stop(
        "Provide either metrics_input or one or more existing single_dir paths.",
        call. = FALSE
      )
    }

    result_files <- unlist(
      lapply(single_dir, function(current_dir) {
        list.files(
          current_dir,
          pattern = "^PHY_kfoldcv_.+_lambda_TTF\\.rds$",
          full.names = TRUE
        )
      }),
      use.names = FALSE
    )

    if (length(result_files) == 0L) {
      stop("No single-trait lambda-TTF PHY files were found.", call. = FALSE)
    }

    # Directory order defines priority. This allows a primary result folder to
    # be supplemented by an older folder without assessing duplicated traits
    # twice (for example, Root depth and Seed dry mass from focal-eight runs).
    result_basenames <- basename(result_files)
    result_files <- result_files[!duplicated(result_basenames)]
    names(result_files) <- tools::file_path_sans_ext(basename(result_files))

    if (!exists("calculate_phy_metrics", mode = "function")) {
      if (!file.exists(calculate_metrics_file)) {
        stop("calculate_phy_metrics.R was not found.", call. = FALSE)
      }
      source(calculate_metrics_file, local = environment())
    }

    metrics <- calculate_phy_metrics(
      phy_results = result_files,
      traits = traits,
      intervals = interval,
      use_transformed_truth = use_transformed_truth,
      folds = list(lambda = folds),
      save_dir = NULL
    )
  }

  if (!is.list(metrics) || !is.data.frame(metrics$predictions)) {
    stop("metrics_input must contain a data-frame element named 'predictions'.", call. = FALSE)
  }

  predictions <- metrics$predictions
  required <- c(
    "result_name", "fold", "species", "trait", "truth", "prediction"
  )
  missing_required <- setdiff(required, names(predictions))
  if (length(missing_required) > 0L) {
    stop(
      paste("The predictions table is missing:", paste(missing_required, collapse = ", ")),
      call. = FALSE
    )
  }

  interval_label <- gsub("\\.", "_", format(interval, trim = TRUE, scientific = FALSE))
  lower_column <- paste0("lower_", interval_label)
  upper_column <- paste0("upper_", interval_label)
  covered_column <- paste0("covered_", interval_label)

  missing_interval <- setdiff(
    c(lower_column, upper_column, covered_column),
    names(predictions)
  )
  if (length(missing_interval) > 0L) {
    stop(
      paste0(
        "The requested interval is absent from metrics_input: ",
        paste(missing_interval, collapse = ", "),
        ". Recalculate metrics with intervals = ", interval, "."
      ),
      call. = FALSE
    )
  }

  is_single_lambda_ttf <- (
    grepl("^PHY_kfoldcv_.+_lambda_TTF$", predictions$result_name) &
      !grepl("all8|all_traits|high_correlation", predictions$result_name)
  )
  predictions <- predictions[is_single_lambda_ttf, , drop = FALSE]

  if (!is.null(traits)) {
    predictions <- predictions[predictions$trait %in% traits, , drop = FALSE]
  }

  if (nrow(predictions) == 0L) {
    stop("No single-trait lambda-TTF predictions were found.", call. = FALSE)
  }

  predictions$lower <- as.numeric(predictions[[lower_column]])
  predictions$upper <- as.numeric(predictions[[upper_column]])
  predictions$covered <- as.logical(predictions[[covered_column]])
  predictions$interval_width <- predictions$upper - predictions$lower

  finite_rows <- with(
    predictions,
    is.finite(truth) & is.finite(prediction) &
      is.finite(lower) & is.finite(upper) & !is.na(covered)
  )
  predictions <- predictions[finite_rows, , drop = FALSE]

  if (anyDuplicated(predictions[c("result_name", "fold", "species", "trait")])) {
    stop(
      paste0(
        "Duplicated fold-species-trait assessment units were found. ",
        "Do not replicate a species mean across its original individual records."
      ),
      call. = FALSE
    )
  }

  trait_order <- if (is.null(traits)) {
    sort(unique(as.character(predictions$trait)))
  } else {
    intersect(traits, unique(as.character(predictions$trait)))
  }

  summary_rows <- lapply(trait_order, function(trait_name) {
    current <- predictions[predictions$trait == trait_name, , drop = FALSE]
    data.frame(
      trait = trait_name,
      n_predictions = nrow(current),
      coverage = 100 * mean(current$covered),
      median_interval_width = stats::median(current$interval_width),
      median_absolute_error = stats::median(abs(current$prediction - current$truth)),
      stringsAsFactors = FALSE
    )
  })
  coverage_summary <- do.call(rbind, summary_rows)

  predictions$trait_label <- factor(
    predictions$trait,
    levels = trait_order,
    labels = gsub("_", " ", trait_order, fixed = TRUE)
  )

  # Deterministic display-only thinning prevents thousands of intervals from
  # becoming an opaque block. Coverage labels always use every prediction.
  set.seed(as.integer(display_seed))
  display_indices <- unlist(lapply(trait_order, function(trait_name) {
    indices <- which(predictions$trait == trait_name)
    if (is.infinite(max_points_per_trait) || length(indices) <= max_points_per_trait) {
      indices
    } else {
      sort(sample(indices, as.integer(max_points_per_trait), replace = FALSE))
    }
  }), use.names = FALSE)
  display_data <- predictions[display_indices, , drop = FALSE]

  pages <- split(
    trait_order,
    ceiling(seq_along(trait_order) / traits_per_page)
  )

  status_colors <- c(
    "Covered" = "#2C7FB8",
    "Not covered" = "#D95F0E"
  )

  plots <- vector("list", length(pages))
  saved_files <- character(0)

  for (page_index in seq_along(pages)) {
    page_traits <- pages[[page_index]]
    panel_plots <- lapply(page_traits, function(trait_name) {
      current <- display_data[display_data$trait == trait_name, , drop = FALSE]
      current$status <- factor(
        ifelse(current$covered, "Covered", "Not covered"),
        levels = c("Covered", "Not covered")
      )

      covered_data <- current[current$covered, , drop = FALSE]
      missed_data <- current[!current$covered, , drop = FALSE]

      limits <- range(
        c(current$truth, current$prediction, current$lower, current$upper),
        finite = TRUE
      )
      padding <- max(diff(limits) * 0.035, 0.05)
      limits <- limits + c(-padding, padding)

      n_current <- nrow(current)
      covered_line_alpha <- max(0.018, min(0.16, 120 / n_current))
      missed_line_alpha <- max(0.055, min(0.30, 220 / n_current))
      covered_point_alpha <- max(0.08, min(0.42, 260 / n_current))
      missed_point_alpha <- max(0.16, min(0.60, 420 / n_current))

      current_summary <- coverage_summary[
        coverage_summary$trait == trait_name,
        ,
        drop = FALSE
      ]
      panel_title <- gsub("_", " ", trait_name, fixed = TRUE)
      if (identical(trait_name, "Leaf_Vcmax_per_dry_mass")) {
        panel_title <- "Leaf Vcmax per\ndry mass"
      }
      panel_subtitle <- sprintf(
        "Coverage %.1f%%\nn = %s",
        current_summary$coverage,
        format(current_summary$n_predictions, big.mark = ",", scientific = FALSE)
      )

      ggplot2::ggplot() +
        ggplot2::geom_abline(
          intercept = 0,
          slope = 1,
          color = "grey18",
          linewidth = 0.38
        ) +
        ggplot2::geom_segment(
          data = covered_data,
          ggplot2::aes(
            x = truth,
            xend = truth,
            y = lower,
            yend = upper,
            color = status
          ),
          linewidth = 0.22,
          alpha = covered_line_alpha,
          lineend = "round"
        ) +
        ggplot2::geom_point(
          data = covered_data,
          ggplot2::aes(x = truth, y = prediction, color = status),
          size = 0.48,
          alpha = covered_point_alpha
        ) +
        ggplot2::geom_segment(
          data = missed_data,
          ggplot2::aes(
            x = truth,
            xend = truth,
            y = lower,
            yend = upper,
            color = status
          ),
          linewidth = 0.24,
          alpha = missed_line_alpha,
          lineend = "round"
        ) +
        ggplot2::geom_point(
          data = missed_data,
          ggplot2::aes(x = truth, y = prediction, color = status),
          size = 0.52,
          alpha = missed_point_alpha
        ) +
        ggplot2::scale_color_manual(values = status_colors, drop = FALSE) +
        ggplot2::scale_x_continuous(limits = limits) +
        ggplot2::scale_y_continuous(limits = limits) +
        # Do not force square panels here. With 18 panels on one manuscript
        # page, coord_equal() shrinks the entire four-column grid and creates
        # large unused margins. Equal x/y limits still make the black line
        # represent y = x, while the panels can now occupy the full width.
        ggplot2::coord_cartesian(expand = FALSE) +
        ggplot2::labs(
          title = panel_title,
          subtitle = panel_subtitle,
          x = NULL,
          y = NULL,
          color = NULL
        ) +
        ggplot2::theme_bw(base_size = base_size) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = base_size + 0.3,
            face = "plain",
            hjust = 0,
            margin = ggplot2::margin(0, 0, 1, 0)
          ),
          plot.subtitle = ggplot2::element_text(
            size = base_size - 0.5,
            color = "grey38",
            hjust = 0,
            lineheight = 0.95,
            margin = ggplot2::margin(0, 0, 5, 0)
          ),
          axis.text = ggplot2::element_text(size = base_size - 0.7, color = "grey25"),
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.22),
          panel.border = ggplot2::element_rect(color = "grey55", linewidth = 0.35),
          legend.position = "none",
          # Keep only a modest inter-column gutter so the plotting rectangles
          # occupy as much of the fixed-width manuscript canvas as possible.
          plot.margin = ggplot2::margin(4, 1, 2, 1)
        )
    })

    # The manuscript figure contains 18 traits on one page. The first 16
    # panels fill a 4 x 4 grid; panels 17 and 18 are centred in the final row.
    # Other page sizes retain the ordinary four-column layout.
    if (length(panel_plots) == 18L) {
      main_panels <- patchwork::wrap_plots(
        panel_plots,
        design = paste(
          "ABCD",
          "EFGH",
          "IJKL",
          "MNOP",
          "#QR#",
          sep = "\n"
        )
      )
    } else {
      main_panels <- patchwork::wrap_plots(
        panel_plots,
        ncol = 4
      )
    }

    legend_panel <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = 0.20, y = 0.5,
        label = "Interval status", hjust = 1,
        size = base_size / ggplot2::.pt,
        color = "grey25"
      ) +
      ggplot2::annotate(
        "segment", x = 0.25, xend = 0.32, y = 0.5, yend = 0.5,
        linewidth = 0.45, color = status_colors[["Covered"]]
      ) +
      ggplot2::annotate(
        "point", x = 0.285, y = 0.5,
        size = 1.4, color = status_colors[["Covered"]]
      ) +
      ggplot2::annotate(
        "text", x = 0.335, y = 0.5,
        label = "Covered", hjust = 0,
        size = base_size / ggplot2::.pt,
        color = "grey25"
      ) +
      ggplot2::annotate(
        "segment", x = 0.52, xend = 0.59, y = 0.5, yend = 0.5,
        linewidth = 0.45, color = status_colors[["Not covered"]]
      ) +
      ggplot2::annotate(
        "point", x = 0.555, y = 0.5,
        size = 1.4, color = status_colors[["Not covered"]]
      ) +
      ggplot2::annotate(
        "text", x = 0.605, y = 0.5,
        label = "Not covered", hjust = 0,
        size = base_size / ggplot2::.pt,
        color = "grey25"
      ) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
      ggplot2::theme_void()

    x_axis_panel <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = "Observed species mean",
        size = base_size / ggplot2::.pt,
        color = "grey20"
      ) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
      ggplot2::theme_void()

    # Draw the shared y-axis title directly on the full panel assembly. It
    # therefore occupies no grid column, remains close to the first column,
    # and leaves the four plotting columns the full available width.
    labelled_panels <- cowplot::ggdraw(main_panels) +
      cowplot::draw_label(
        "Predicted species mean",
        # The aligned patchwork reserves space for the tick labels of the
        # first column. Position the shared label immediately outside that
        # reserved strip instead of at the edge of the exported canvas.
        x = 0.018,
        y = 0.50,
        angle = 90,
        hjust = 0.5,
        vjust = 0.5,
        size = base_size,
        color = "grey20"
      )

    x_axis_row <- x_axis_panel

    legend_row <- legend_panel

    plot_core <-
      labelled_panels /
      x_axis_row /
      legend_row +
      patchwork::plot_layout(heights = c(1, 0.024, 0.035))

    # Preserve the same internal left and right whitespace as the original
    # 4 x 4 manuscript figure. Padding is part of the exported image, so the
    # figure remains centred when Word displays it at the full text width.
    plot_object <-
      (
        patchwork::plot_spacer() |
          plot_core |
          patchwork::plot_spacer()
      ) +
      patchwork::plot_layout(
        widths = c(
          horizontal_padding,
          1 - 2 * horizontal_padding,
          horizontal_padding
        )
      )

    plots[[page_index]] <- plot_object

    if (!is.null(output_dir)) {
      output_file <- file.path(
        output_dir,
        sprintf("%s_page_%02d.png", file_prefix, page_index)
      )
      ggplot2::ggsave(
        filename = output_file,
        plot = plot_object,
        width = width,
        height = height,
        units = "in",
        dpi = dpi,
        bg = "white"
      )
      saved_files <- c(saved_files, output_file)
    }
  }

  if (!is.null(output_dir)) {
    summary_file <- file.path(output_dir, paste0(file_prefix, "_summary.csv"))
    utils::write.csv(coverage_summary, summary_file, row.names = FALSE)
    saved_files <- c(saved_files, summary_file)
  }

  if (interactive()) {
    print(plots[[1L]])
  }

  invisible(list(
    plots = plots,
    coverage_summary = coverage_summary,
    predictions = predictions,
    display_data = display_data,
    saved_files = saved_files
  ))
}


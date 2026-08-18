calculate_phy_metrics <- function(
    phy_results,
    traits = NULL,
    intervals = c(90, 95, 99.5),
    use_transformed_truth = TRUE,
    folds = NULL,
    save_dir = NULL,
    save_name = "phy_metrics"
) {
  phy_expression <- deparse(
    substitute(phy_results),
    width.cutoff = 500L
  )

  if (
    length(use_transformed_truth) != 1L ||
    !is.logical(use_transformed_truth) ||
    is.na(use_transformed_truth)
  ) {
    stop(
      "use_transformed_truth must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  truth_column <- if (use_transformed_truth) {
    "transformed_truth"
  } else {
    "raw_truth"
  }

  if (!is.null(save_dir)) {
    if (
      !is.character(save_dir) ||
      length(save_dir) != 1L ||
      is.na(save_dir) ||
      !nzchar(save_dir)
    ) {
      stop(
        "save_dir must be NULL or one non-empty directory path.",
        call. = FALSE
      )
    }
  }

  if (
    !is.character(save_name) ||
    length(save_name) != 1L ||
    is.na(save_name) ||
    !nzchar(save_name) ||
    basename(save_name) != save_name ||
    save_name %in% c(".", "..")
  ) {
    stop(
      "save_name must be one non-empty file-name prefix.",
      call. = FALSE
    )
  }

  if (!is.null(traits)) {
    if (
      !is.character(traits) ||
      length(traits) == 0L ||
      anyNA(traits) ||
      any(!nzchar(traits))
    ) {
      stop(
        "traits must be NULL or a non-empty character vector.",
        call. = FALSE
      )
    }

    traits <- unique(traits)
  }

  if (
    !is.numeric(intervals) ||
    length(intervals) == 0L ||
    anyNA(intervals) ||
    any(!is.finite(intervals)) ||
    anyDuplicated(intervals)
  ) {
    stop(
      "intervals must contain unique numeric values.",
      call. = FALSE
    )
  }

  if (all(intervals > 1 & intervals < 100)) {
    interval_levels <- intervals / 100
  } else if (all(intervals > 0 & intervals < 1)) {
    interval_levels <- intervals
  } else {
    stop(
      paste0(
        "intervals must be probabilities between 0 and 1 ",
        "or percentages between 1 and 100."
      ),
      call. = FALSE
    )
  }

  interval_levels <- sort(interval_levels)

  normalize_folds <- function(fold_values, argument_name) {
    if (
      !(is.numeric(fold_values) || is.character(fold_values)) ||
      length(fold_values) == 0L ||
      anyNA(fold_values)
    ) {
      stop(
        paste0(
          argument_name,
          " must contain numeric fold numbers or fold names."
        ),
        call. = FALSE
      )
    }

    if (is.numeric(fold_values)) {
      if (
        any(!is.finite(fold_values)) ||
        any(fold_values < 1) ||
        any(fold_values != as.integer(fold_values))
      ) {
        stop(
          paste0(
            argument_name,
            " must contain positive integers."
          ),
          call. = FALSE
        )
      }

      fold_values <- paste0(
        "fold_",
        as.integer(fold_values)
      )
    }

    unique(as.character(fold_values))
  }

  if (!is.null(folds)) {
    if (is.list(folds)) {
      if (
        length(folds) == 0L ||
        is.null(names(folds)) ||
        anyNA(names(folds)) ||
        any(!nzchar(names(folds))) ||
        anyDuplicated(names(folds))
      ) {
        stop(
          paste0(
            "A model-specific folds list must have unique, ",
            "non-empty model names."
          ),
          call. = FALSE
        )
      }

      folds <- stats::setNames(
        lapply(
          seq_along(folds),
          function(fold_index) {
            normalize_folds(
              folds[[fold_index]],
              paste0(
                "folds$",
                names(folds)[[fold_index]]
              )
            )
          }
        ),
        names(folds)
      )
    } else {
      folds <- normalize_folds(
        folds,
        "folds"
      )
    }
  }

  is_phy_result <- inherits(
    phy_results,
    "kfold_phylopars_results"
  )

  if (is_phy_result) {
    result_name <- if (
      length(phy_expression) == 1L &&
      grepl(
        "^[A-Za-z.][A-Za-z0-9._]*$",
        phy_expression
      )
    ) {
      phy_expression
    } else {
      "result_1"
    }

    phy_results <- stats::setNames(
      list(phy_results),
      result_name
    )
  } else if (
    is.character(phy_results) &&
    length(phy_results) == 1L
  ) {
    phy_results <- list(phy_results)
  } else if (
    is.character(phy_results) &&
    length(phy_results) > 1L
  ) {
    phy_results <- as.list(phy_results)
  } else if (!is.list(phy_results)) {
    stop(
      paste0(
        "phy_results must be a PHY result, an RDS path, ",
        "or a list of PHY results or RDS paths."
      ),
      call. = FALSE
    )
  }

  if (length(phy_results) == 0L) {
    stop(
      "phy_results must not be empty.",
      call. = FALSE
    )
  }

  result_names <- names(phy_results)

  if (is.null(result_names)) {
    result_names <- rep(
      "",
      length(phy_results)
    )
  }

  missing_result_names <- !nzchar(result_names)

  if (any(missing_result_names)) {
    result_names[missing_result_names] <- paste0(
      "result_",
      which(missing_result_names)
    )
  }

  if (anyDuplicated(result_names)) {
    stop(
      "phy_results must have unique names.",
      call. = FALSE
    )
  }

  level_labels <- vapply(
    interval_levels,
    function(level) {
      label <- format(
        level * 100,
        trim = TRUE,
        scientific = FALSE,
        digits = 6
      )

      gsub(
        ".",
        "_",
        label,
        fixed = TRUE
      )
    },
    character(1)
  )

  z_scores <- stats::qnorm(
    (1 + interval_levels) / 2
  )

  prediction_tables <- vector(
    "list",
    length(phy_results)
  )

  result_information <- vector(
    "list",
    length(phy_results)
  )

  for (result_index in seq_along(phy_results)) {
    result_input <- phy_results[[result_index]]
    supplied_name <- result_names[[result_index]]

    if (
      is.character(result_input) &&
      length(result_input) == 1L
    ) {
      if (!file.exists(result_input)) {
        stop(
          paste0(
            "PHY result file does not exist: ",
            result_input
          ),
          call. = FALSE
        )
      }

      result <- tryCatch(
        readRDS(result_input),
        error = function(error) error
      )

      if (inherits(result, "error")) {
        stop(
          paste0(
            "Could not read PHY result file: ",
            result_input
          ),
          call. = FALSE
        )
      }

      if (grepl("^result_[0-9]+$", supplied_name)) {
        supplied_name <- tools::file_path_sans_ext(
          basename(result_input)
        )
      }
    } else {
      result <- result_input
    }

    message(
      sprintf(
        "[%d/%d] Calculating: %s",
        result_index,
        length(phy_results),
        supplied_name
      )
    )

    status_table <- attr(result, "status")
    assessments <- attr(result, "assessment")

    valid_result <- (
      inherits(
        result,
        "kfold_phylopars_results"
      ) &&
        is.list(result) &&
        !is.null(names(result)) &&
        is.data.frame(status_table) &&
        all(
          c("fold", "status") %in%
            names(status_table)
        ) &&
        is.list(assessments) &&
        length(assessments) == length(result) &&
        nrow(status_table) == length(result)
    )

    if (!valid_result) {
      stop(
        paste0(
          "Invalid PHY result: ",
          supplied_name
        ),
        call. = FALSE
      )
    }

    available_fold_names <- as.character(
      status_table$fold
    )

    result_model <- attr(result, "model")

    if (
      is.null(result_model) ||
      length(result_model) == 0L
    ) {
      result_model <- as.character(
        status_table$model[[1L]]
      )
    } else {
      result_model <- as.character(
        result_model[[1L]]
      )
    }

    selected_fold_names <- if (is.null(folds)) {
      available_fold_names
    } else if (is.list(folds)) {
      if (!(result_model %in% names(folds))) {
        stop(
          paste0(
            "No fold selection was supplied for model: ",
            result_model
          ),
          call. = FALSE
        )
      }

      folds[[result_model]]
    } else {
      folds
    }

    missing_folds <- setdiff(
      selected_fold_names,
      available_fold_names
    )

    if (length(missing_folds) > 0L) {
      stop(
        paste0(
          "Requested folds not found in ",
          supplied_name,
          ": ",
          paste(missing_folds, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    selected_fold_indices <- match(
      selected_fold_names,
      available_fold_names
    )

    config_match <- regexec(
      "_(BM|lambda)_([FT]{3})$",
      supplied_name
    )

    config_parts <- regmatches(
      supplied_name,
      config_match
    )[[1L]]

    configuration <- if (length(config_parts) >= 3L) {
      config_parts[[3L]]
    } else {
      NA_character_
    }

    result_information[[result_index]] <- data.frame(
      result_name = supplied_name,
      model = result_model,
      configuration = configuration,
      stringsAsFactors = FALSE
    )

    fold_tables <- vector(
      "list",
      length(selected_fold_indices)
    )

    for (
      selected_index in
      seq_along(selected_fold_indices)
    ) {
      fold_index <- selected_fold_indices[[selected_index]]

      if (
        !identical(
          as.character(
            status_table$status[[fold_index]]
          ),
          "success"
        ) ||
        is.null(result[[fold_index]])
      ) {
        next
      }

      fold_result <- result[[fold_index]]
      assessment <- assessments[[fold_index]]

      valid_fold <- (
        is.list(fold_result) &&
          all(
            c("anc_recon", "anc_var") %in%
              names(fold_result)
          ) &&
          is.data.frame(assessment) &&
          all(
            c(
              "species",
              "trait",
              truth_column
            ) %in% names(assessment)
          )
      )

      if (!valid_fold) {
        stop(
          paste0(
            "Invalid fold data in: ",
            supplied_name
          ),
          call. = FALSE
        )
      }

      available_traits <- intersect(
        colnames(fold_result$anc_recon),
        colnames(fold_result$anc_var)
      )

      if (!is.null(traits)) {
        available_traits <- intersect(
          available_traits,
          traits
        )
      }

      assessment <- assessment[
        is.finite(assessment[[truth_column]]) &
          assessment$trait %in% available_traits,
        ,
        drop = FALSE
      ]

      if (nrow(assessment) == 0L) {
        next
      }

      assessment_groups <- data.frame(
        species = as.character(assessment$species),
        trait = as.character(assessment$trait),
        truth = assessment[[truth_column]],
        stringsAsFactors = FALSE
      )

      grouped_truth <- stats::aggregate(
        truth ~ species + trait,
        data = assessment_groups,
        FUN = mean
      )

      recon_rows <- match(
        grouped_truth$species,
        rownames(fold_result$anc_recon)
      )

      recon_columns <- match(
        grouped_truth$trait,
        colnames(fold_result$anc_recon)
      )

      variance_rows <- match(
        grouped_truth$species,
        rownames(fold_result$anc_var)
      )

      variance_columns <- match(
        grouped_truth$trait,
        colnames(fold_result$anc_var)
      )

      if (
        anyNA(recon_rows) ||
        anyNA(recon_columns) ||
        anyNA(variance_rows) ||
        anyNA(variance_columns)
      ) {
        stop(
          paste0(
            "Prediction indices are incomplete in: ",
            supplied_name
          ),
          call. = FALSE
        )
      }

      prediction <- as.numeric(
        fold_result$anc_recon[
          cbind(recon_rows, recon_columns)
        ]
      )

      prediction_variance <- as.numeric(
        fold_result$anc_var[
          cbind(variance_rows, variance_columns)
        ]
      )

      if (
        any(!is.finite(prediction)) ||
        any(!is.finite(prediction_variance))
      ) {
        stop(
          paste0(
            "Non-finite predictions found in: ",
            supplied_name
          ),
          call. = FALSE
        )
      }

      if (
        any(
          prediction_variance <
          -sqrt(.Machine$double.eps)
        )
      ) {
        stop(
          paste0(
            "Negative prediction variance found in: ",
            supplied_name
          ),
          call. = FALSE
        )
      }

      prediction_variance <- pmax(
        prediction_variance,
        0
      )

      prediction_se <- sqrt(
        prediction_variance
      )

      truth <- grouped_truth$truth

      output <- data.frame(
        result_name = supplied_name,
        model = result_model,
        configuration = configuration,
        fold = available_fold_names[[fold_index]],
        species = grouped_truth$species,
        trait = grouped_truth$trait,
        truth = truth,
        prediction = prediction,
        prediction_variance = prediction_variance,
        prediction_se = prediction_se,
        residual = prediction - truth,
        absolute_error = abs(prediction - truth),
        squared_error = (prediction - truth)^2,
        stringsAsFactors = FALSE
      )

      for (level_index in seq_along(interval_levels)) {
        label <- level_labels[[level_index]]
        margin <- z_scores[[level_index]] * prediction_se
        lower <- prediction - margin
        upper <- prediction + margin

        output[[paste0("lower_", label)]] <- lower
        output[[paste0("upper_", label)]] <- upper
        output[[paste0("covered_", label)]] <- (
          truth >= lower &
            truth <= upper
        )
      }

      fold_tables[[selected_index]] <- output
    }

    fold_tables <- fold_tables[
      !vapply(
        fold_tables,
        is.null,
        logical(1)
      )
    ]

    if (length(fold_tables) > 0L) {
      prediction_tables[[result_index]] <- do.call(
        rbind,
        fold_tables
      )
    }

    message(
      sprintf(
        "[%d/%d] Calculation completed: %s",
        result_index,
        length(phy_results),
        supplied_name
      )
    )

    rm(result)
    invisible(gc(verbose = FALSE))
  }

  prediction_tables <- prediction_tables[
    !vapply(
      prediction_tables,
      is.null,
      logical(1)
    )
  ]

  if (length(prediction_tables) == 0L) {
    stop(
      "No successful predictions were found.",
      call. = FALSE
    )
  }

  predictions <- do.call(
    rbind,
    prediction_tables
  )

  rownames(predictions) <- NULL

  if (!is.null(traits)) {
    missing_traits <- setdiff(
      traits,
      unique(predictions$trait)
    )

    if (length(missing_traits) > 0L) {
      stop(
        paste0(
          "Requested traits were not found: ",
          paste(missing_traits, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  metric_key <- interaction(
    predictions$result_name,
    predictions$trait,
    drop = TRUE,
    lex.order = TRUE
  )

  metric_indices <- split(
    seq_len(nrow(predictions)),
    metric_key
  )

  metric_rows <- lapply(
    metric_indices,
    function(indices) {
      current <- predictions[
        indices,
        ,
        drop = FALSE
      ]

      truth <- current$truth
      prediction <- current$prediction
      residual <- prediction - truth
      denominator <- sum(
        (truth - mean(truth))^2
      )

      output <- data.frame(
        result_name = current$result_name[[1L]],
        model = current$model[[1L]],
        configuration = current$configuration[[1L]],
        trait = current$trait[[1L]],
        folds_used = paste(
          unique(current$fold),
          collapse = ", "
        ),
        rmse = sqrt(mean(residual^2)),
        mae = mean(abs(residual)),
        r2 = if (denominator > 0) {
          1 - sum(residual^2) / denominator
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )

      for (level_index in seq_along(interval_levels)) {
        label <- level_labels[[level_index]]

        output[[paste0("coverage_", label)]] <- (
          100 * mean(
            current[[paste0("covered_", label)]]
          )
        )
      }

      output
    }
  )

  trait_metrics <- do.call(
    rbind,
    metric_rows
  )

  rownames(trait_metrics) <- NULL

  trait_metrics <- trait_metrics[
    order(
      trait_metrics$result_name,
      trait_metrics$trait
    ),
    ,
    drop = FALSE
  ]

  predictions <- predictions[
    order(
      predictions$result_name,
      predictions$trait,
      predictions$fold,
      predictions$species
    ),
    ,
    drop = FALSE
  ]

  rownames(trait_metrics) <- NULL
  rownames(predictions) <- NULL

  metric_columns <- c(
    "rmse",
    "mae",
    "r2",
    grep(
      "^coverage_",
      names(trait_metrics),
      value = TRUE
    )
  )

  result_information <- do.call(
    rbind,
    result_information
  )

  rownames(result_information) <- NULL

  result_information$analysis_type <- vapply(
    result_information$result_name,
    function(result_name) {
      if (grepl("(^|_)all_traits_", result_name)) {
        "all_traits"
      } else if (
        grepl(
          "(^|_)high_correlation_traits_",
          result_name
        )
      ) {
        "high_correlation_traits"
      } else if (
        grepl(
          "(^|_)high_correlation_single_",
          result_name
        )
      ) {
        "high_correlation_single"
      } else if (grepl("(^|_)all8_", result_name)) {
        "all8"
      } else {
        "single"
      }
    },
    character(1)
  )

  result_information$analysis <- paste(
    result_information$analysis_type,
    result_information$model,
    result_information$configuration,
    sep = "_"
  )

  metric_source <- merge(
    trait_metrics,
    result_information[
      c(
        "result_name",
        "analysis_type",
        "analysis"
      )
    ],
    by = "result_name",
    all.x = TRUE,
    sort = FALSE
  )

  duplicate_metric_keys <- paste(
    metric_source$analysis,
    metric_source$trait,
    sep = "\r"
  )

  if (anyDuplicated(duplicate_metric_keys)) {
    stop(
      paste0(
        "Multiple PHY results supply the same analysis and trait. ",
        "Use unique result names and configurations."
      ),
      call. = FALSE
    )
  }

  single_trait_columns <- sort(
    unique(
      as.character(
        metric_source$trait[
          metric_source$analysis_type == "single"
        ]
      )
    )
  )

  all_trait_columns <- sort(
    unique(
      as.character(metric_source$trait)
    )
  )

  trait_columns <- c(
    single_trait_columns,
    setdiff(
      all_trait_columns,
      single_trait_columns
    )
  )

  analysis_information <- unique(
    result_information[
      c(
        "analysis",
        "analysis_type",
        "model",
        "configuration"
      )
    ]
  )

  analysis_information <- analysis_information[
    order(
      match(
        analysis_information$analysis_type,
        c(
          "single",
          "all8",
          "high_correlation_single",
          "high_correlation_traits",
          "all_traits"
        )
      ),
      match(
        analysis_information$model,
        c("BM", "lambda")
      ),
      match(
        analysis_information$configuration,
        c("FTF", "TTF")
      ),
      analysis_information$analysis
    ),
    ,
    drop = FALSE
  ]

  rownames(analysis_information) <- NULL

  metric_tables <- stats::setNames(
    vector(
      "list",
      length(metric_columns)
    ),
    metric_columns
  )

  for (metric_column in metric_columns) {
    metric_table <- data.frame(
      trait = trait_columns,
      stringsAsFactors = FALSE
    )

    for (analysis_column in analysis_information$analysis) {
      metric_table[[analysis_column]] <- NA_real_
    }

    for (
      analysis_information_index in seq_len(
        nrow(analysis_information)
      )
    ) {
      current_analysis <- analysis_information$analysis[[
        analysis_information_index
      ]]

      current_metrics <- metric_source[
        metric_source$analysis == current_analysis,
        ,
        drop = FALSE
      ]

      if (nrow(current_metrics) == 0L) {
        next
      }

      trait_positions <- match(
        current_metrics$trait,
        trait_columns
      )

      metric_table[[current_analysis]][
        trait_positions
      ] <- current_metrics[[metric_column]]
    }

    rownames(metric_table) <- NULL
    metric_tables[[metric_column]] <- metric_table
  }

  results <- structure(
    list(
      trait_metrics = trait_metrics,
      metric_tables = metric_tables,
      analysis_information = analysis_information,
      predictions = predictions,
      intervals = interval_levels,
      z_scores = stats::setNames(
        z_scores,
        level_labels
      ),
      truth_column = truth_column,
      selected_folds = folds,
      selected_traits = traits
    ),
    class = "phy_metric_results"
  )

  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) {
      directory_created <- dir.create(
        save_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )

      if (!directory_created && !dir.exists(save_dir)) {
        stop(
          paste0(
            "Could not create save directory: ",
            save_dir
          ),
          call. = FALSE
        )
      }
    }

    rds_path <- file.path(
      save_dir,
      paste0(save_name, ".rds")
    )

    metric_csv_paths <- stats::setNames(
      file.path(
        save_dir,
        paste0(
          save_name,
          "_",
          metric_columns,
          ".csv"
        )
      ),
      metric_columns
    )

    results$saved_files <- c(
      rds = rds_path,
      stats::setNames(
        metric_csv_paths,
        paste0(
          names(metric_csv_paths),
          "_csv"
        )
      )
    )

    saveRDS(
      results,
      rds_path
    )

    for (metric_column in metric_columns) {
      utils::write.csv(
        metric_tables[[metric_column]],
        metric_csv_paths[[metric_column]],
        row.names = FALSE,
        na = ""
      )

      message(
        sprintf(
          "Metric table saved (%s): %s",
          metric_column,
          metric_csv_paths[[metric_column]]
        )
      )
    }

    message(
      paste0(
        "Results saved: ",
        rds_path
      )
    )
  }

  results
}

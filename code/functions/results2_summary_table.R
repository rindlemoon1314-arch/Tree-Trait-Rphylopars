build_results2_summary_table <- function(
    trait_data,
    tree,
    r2_data,
    trait_groups,
    group_labels = NULL,
    species_column = "species",
    configurations = c(
      "BM_FTF",
      "BM_TTF",
      "lambda_FTF",
      "lambda_TTF"
    ),
    single_label = "Single trait",
    multi_label = "Multi-trait",
    delta_tolerance = 1e-12,
    output_dir = NULL,
    file_prefix = "Results2_single_multi",
    save_xlsx = TRUE
) {
  stop_with <- function(...) {
    stop(sprintf(...), call. = FALSE)
  }

  require_columns <- function(data, columns, object_name) {
    missing_columns <- setdiff(columns, names(data))
    if (length(missing_columns) > 0L) {
      stop_with(
        "%s is missing required columns: %s",
        object_name,
        paste(missing_columns, collapse = ", ")
      )
    }
  }

  read_table_input <- function(x, object_name) {
    if (is.data.frame(x)) {
      return(x)
    }

    if (is.character(x) && length(x) == 1L && file.exists(x)) {
      return(
        utils::read.csv(
          x,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      )
    }

    stop_with(
      "%s must be a data frame or the path to a readable CSV file.",
      object_name
    )
  }

  extract_tip_labels <- function(x) {
    if (inherits(x, "phylo")) {
      return(as.character(x$tip.label))
    }

    if (is.character(x) && length(x) == 1L && file.exists(x)) {
      if (!requireNamespace("ape", quietly = TRUE)) {
        stop_with(
          "Package 'ape' is required when tree is supplied as a file path."
        )
      }

      return(as.character(ape::read.tree(x)$tip.label))
    }

    if (is.character(x) && length(x) >= 1L) {
      return(as.character(x))
    }

    stop_with(
      paste0(
        "tree must be a phylo object, a Newick file path, ",
        "or a character vector of tip labels."
      )
    )
  }

  humanize_trait <- function(x) {
    gsub("_", " ", x, fixed = TRUE)
  }

  format_integer_range <- function(minimum, maximum) {
    if (is.na(minimum) || is.na(maximum)) {
      return(NA_character_)
    }

    if (minimum == maximum) {
      return(format(minimum, scientific = FALSE, trim = TRUE))
    }

    paste0(
      format(minimum, scientific = FALSE, trim = TRUE),
      "-",
      format(maximum, scientific = FALSE, trim = TRUE)
    )
  }

  write_csv_file <- function(data, path) {
    utils::write.csv(
      data,
      path,
      row.names = FALSE,
      na = ""
    )
  }

  write_excel_file <- function(tables, path) {
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      workbook <- openxlsx::createWorkbook()

      header_style <- openxlsx::createStyle(
        fgFill = "#1F4E78",
        fontColour = "#FFFFFF",
        textDecoration = "bold",
        halign = "center",
        valign = "center",
        wrapText = TRUE,
        border = "Bottom",
        borderColour = "#A6A6A6"
      )

      for (sheet_name in names(tables)) {
        data <- tables[[sheet_name]]
        openxlsx::addWorksheet(workbook, sheet_name)
        openxlsx::writeData(
          workbook,
          sheet = sheet_name,
          x = data,
          withFilter = TRUE
        )
        openxlsx::addStyle(
          workbook,
          sheet = sheet_name,
          style = header_style,
          rows = 1L,
          cols = seq_len(ncol(data)),
          gridExpand = TRUE
        )
        openxlsx::freezePane(
          workbook,
          sheet = sheet_name,
          firstActiveRow = 2L
        )
        openxlsx::setColWidths(
          workbook,
          sheet = sheet_name,
          cols = seq_len(ncol(data)),
          widths = "auto"
        )
      }

      openxlsx::saveWorkbook(
        workbook,
        file = path,
        overwrite = TRUE
      )
      return(TRUE)
    }

    if (requireNamespace("writexl", quietly = TRUE)) {
      writexl::write_xlsx(tables, path)
      return(TRUE)
    }

    warning(
      paste0(
        "The Excel workbook was not written because neither 'openxlsx' ",
        "nor 'writexl' is installed. CSV and RDS files were still saved."
      ),
      call. = FALSE
    )
    FALSE
  }

  trait_data <- read_table_input(trait_data, "trait_data")
  r2_data <- read_table_input(r2_data, "r2_data")

  if (!is.list(trait_groups) || is.null(names(trait_groups))) {
    stop_with("trait_groups must be a named list of character vectors.")
  }

  if (any(names(trait_groups) == "") || anyDuplicated(names(trait_groups))) {
    stop_with("trait_groups must have unique, non-empty names.")
  }

  invalid_groups <- vapply(
    trait_groups,
    function(x) !is.character(x) || length(x) == 0L || anyNA(x),
    logical(1)
  )

  if (any(invalid_groups)) {
    stop_with(
      "Every element of trait_groups must be a non-empty character vector."
    )
  }

  if (is.null(group_labels)) {
    group_labels <- stats::setNames(names(trait_groups), names(trait_groups))
  } else {
    if (is.null(names(group_labels))) {
      stop_with("group_labels must be a named character vector.")
    }

    missing_labels <- setdiff(names(trait_groups), names(group_labels))
    if (length(missing_labels) > 0L) {
      stop_with(
        "group_labels is missing labels for: %s",
        paste(missing_labels, collapse = ", ")
      )
    }

    group_labels <- group_labels[names(trait_groups)]
  }

  all_traits <- unique(unlist(trait_groups, use.names = FALSE))
  require_columns(
    trait_data,
    c(species_column, all_traits),
    "trait_data"
  )

  require_columns(
    r2_data,
    c(
      "group",
      "trait",
      "analysis",
      "configuration_label",
      "folds_used",
      "r2"
    ),
    "r2_data"
  )

  tip_labels <- unique(extract_tip_labels(tree))
  tip_labels <- tip_labels[!is.na(tip_labels) & nzchar(tip_labels)]

  if (length(tip_labels) == 0L) {
    stop_with("No valid tip labels were found in tree.")
  }

  species_values <- as.character(trait_data[[species_column]])
  keep_rows <- !is.na(species_values) &
    nzchar(species_values) &
    species_values %in% tip_labels

  filtered_traits <- trait_data[keep_rows, , drop = FALSE]
  filtered_traits[[species_column]] <- as.character(
    filtered_traits[[species_column]]
  )

  if (nrow(filtered_traits) == 0L) {
    stop_with("No trait_data rows matched the tree tip labels.")
  }

  r2_data$group <- as.character(r2_data$group)
  r2_data$trait <- as.character(r2_data$trait)
  r2_data$analysis <- as.character(r2_data$analysis)
  r2_data$configuration_label <- as.character(
    r2_data$configuration_label
  )
  r2_data$folds_used <- as.character(r2_data$folds_used)
  r2_data$r2 <- suppressWarnings(as.numeric(r2_data$r2))

  if (anyNA(r2_data$r2)) {
    stop_with("r2_data$r2 contains missing or non-numeric values.")
  }

  unknown_groups <- setdiff(unique(r2_data$group), names(trait_groups))
  if (length(unknown_groups) > 0L) {
    stop_with(
      "r2_data contains groups not listed in trait_groups: %s",
      paste(unknown_groups, collapse = ", ")
    )
  }

  missing_r2_groups <- setdiff(names(trait_groups), unique(r2_data$group))
  if (length(missing_r2_groups) > 0L) {
    stop_with(
      "r2_data is missing groups: %s",
      paste(missing_r2_groups, collapse = ", ")
    )
  }

  duplicate_key <- paste(
    r2_data$group,
    r2_data$trait,
    r2_data$analysis,
    r2_data$configuration_label,
    sep = "\r"
  )

  if (anyDuplicated(duplicate_key)) {
    duplicated_rows <- unique(duplicate_key[duplicated(duplicate_key)])
    stop_with(
      paste0(
        "r2_data contains duplicated group-trait-analysis-configuration ",
        "records. First duplicate key: %s"
      ),
      duplicated_rows[1L]
    )
  }

  trait_detail_rows <- list()
  pairwise_rows <- list()
  delta_rows <- list()
  main_rows <- list()

  trait_detail_index <- 0L
  pairwise_index <- 0L
  delta_index <- 0L
  main_index <- 0L

  for (group_key in names(trait_groups)) {
    group_traits <- unique(trait_groups[[group_key]])
    group_label <- unname(group_labels[[group_key]])
    observed_species <- vector("list", length(group_traits))
    names(observed_species) <- group_traits

    trait_observations <- integer(length(group_traits))
    trait_species_counts <- integer(length(group_traits))

    for (i in seq_along(group_traits)) {
      trait_name <- group_traits[i]
      observed <- !is.na(filtered_traits[[trait_name]])
      species_set <- sort(unique(
        filtered_traits[[species_column]][observed]
      ))

      observed_species[[trait_name]] <- species_set
      trait_observations[i] <- sum(observed)
      trait_species_counts[i] <- length(species_set)

      trait_detail_index <- trait_detail_index + 1L
      trait_detail_rows[[trait_detail_index]] <- data.frame(
        group_key = group_key,
        analysis_set = group_label,
        trait_code = trait_name,
        trait = humanize_trait(trait_name),
        observation_records = trait_observations[i],
        observed_species = trait_species_counts[i],
        stringsAsFactors = FALSE
      )
    }

    union_species <- Reduce(union, observed_species)
    shared_species <- Reduce(intersect, observed_species)

    pairwise_shared <- numeric(0L)
    if (length(group_traits) >= 2L) {
      pairs <- utils::combn(group_traits, 2L, simplify = FALSE)

      for (trait_pair in pairs) {
        species_1 <- observed_species[[trait_pair[1L]]]
        species_2 <- observed_species[[trait_pair[2L]]]
        shared_pair <- intersect(species_1, species_2)
        union_pair <- union(species_1, species_2)
        denominator <- min(length(species_1), length(species_2))

        pairwise_shared <- c(pairwise_shared, length(shared_pair))
        pairwise_index <- pairwise_index + 1L
        pairwise_rows[[pairwise_index]] <- data.frame(
          group_key = group_key,
          analysis_set = group_label,
          trait_1_code = trait_pair[1L],
          trait_1 = humanize_trait(trait_pair[1L]),
          trait_2_code = trait_pair[2L],
          trait_2 = humanize_trait(trait_pair[2L]),
          species_trait_1 = length(species_1),
          species_trait_2 = length(species_2),
          shared_species = length(shared_pair),
          union_species = length(union_pair),
          jaccard_index = if (
            length(union_pair) == 0L
          ) NA_real_ else length(shared_pair) / length(union_pair),
          overlap_coefficient = if (
            denominator == 0L
          ) NA_real_ else length(shared_pair) / denominator,
          stringsAsFactors = FALSE
        )
      }
    }

    group_r2 <- r2_data[r2_data$group == group_key, , drop = FALSE]

    unexpected_traits <- setdiff(unique(group_r2$trait), group_traits)
    if (length(unexpected_traits) > 0L) {
      stop_with(
        "r2_data contains unexpected traits for %s: %s",
        group_key,
        paste(unexpected_traits, collapse = ", ")
      )
    }

    missing_traits <- setdiff(group_traits, unique(group_r2$trait))
    if (length(missing_traits) > 0L) {
      stop_with(
        "r2_data is missing traits for %s: %s",
        group_key,
        paste(missing_traits, collapse = ", ")
      )
    }

    for (trait_name in group_traits) {
      for (configuration in configurations) {
        block <- group_r2[
          group_r2$trait == trait_name &
            group_r2$configuration_label == configuration,
          ,
          drop = FALSE
        ]

        single_row <- block[block$analysis == single_label, , drop = FALSE]
        multi_row <- block[block$analysis == multi_label, , drop = FALSE]

        if (nrow(single_row) != 1L || nrow(multi_row) != 1L) {
          stop_with(
            paste0(
              "Expected one '%s' and one '%s' record for ",
              "%s / %s / %s."
            ),
            single_label,
            multi_label,
            group_key,
            trait_name,
            configuration
          )
        }

        if (!identical(single_row$folds_used, multi_row$folds_used)) {
          stop_with(
            "Fold definitions differ between single and multi models for %s / %s / %s.",
            group_key,
            trait_name,
            configuration
          )
        }

        delta_value <- multi_row$r2 - single_row$r2
        direction <- if (
          delta_value > delta_tolerance
        ) {
          "Positive"
        } else if (delta_value < -delta_tolerance) {
          "Negative"
        } else {
          "No change"
        }

        delta_index <- delta_index + 1L
        delta_rows[[delta_index]] <- data.frame(
          group_key = group_key,
          analysis_set = group_label,
          trait_code = trait_name,
          trait = humanize_trait(trait_name),
          configuration = configuration,
          folds_used = single_row$folds_used,
          single_r2 = single_row$r2,
          multi_r2 = multi_row$r2,
          delta_r2 = delta_value,
          direction = direction,
          stringsAsFactors = FALSE
        )
      }
    }

    current_delta <- do.call(
      rbind,
      delta_rows[
        vapply(
          delta_rows,
          function(x) identical(x$group_key[1L], group_key),
          logical(1)
        )
      ]
    )

    pairwise_min <- if (
      length(pairwise_shared) == 0L
    ) NA_integer_ else min(pairwise_shared)
    pairwise_max <- if (
      length(pairwise_shared) == 0L
    ) NA_integer_ else max(pairwise_shared)
    pairwise_mean <- if (
      length(pairwise_shared) == 0L
    ) NA_real_ else mean(pairwise_shared)

    main_index <- main_index + 1L
    main_rows[[main_index]] <- data.frame(
      group_key = group_key,
      analysis_set = group_label,
      traits = paste(humanize_trait(group_traits), collapse = "; "),
      number_of_traits = length(group_traits),
      total_observation_records = sum(trait_observations),
      observations_per_trait_min = min(trait_observations),
      observations_per_trait_max = max(trait_observations),
      observations_per_trait_range = format_integer_range(
        min(trait_observations),
        max(trait_observations)
      ),
      union_species = length(union_species),
      species_per_trait_min = min(trait_species_counts),
      species_per_trait_max = max(trait_species_counts),
      species_per_trait_range = format_integer_range(
        min(trait_species_counts),
        max(trait_species_counts)
      ),
      species_shared_by_all_traits = length(shared_species),
      shared_fraction_of_union = if (
        length(union_species) == 0L
      ) NA_real_ else length(shared_species) / length(union_species),
      pairwise_shared_species_min = pairwise_min,
      pairwise_shared_species_max = pairwise_max,
      pairwise_shared_species_range = format_integer_range(
        pairwise_min,
        pairwise_max
      ),
      pairwise_shared_species_mean = pairwise_mean,
      positive_delta_r2 = sum(current_delta$direction == "Positive"),
      negative_delta_r2 = sum(current_delta$direction == "Negative"),
      no_change_delta_r2 = sum(current_delta$direction == "No change"),
      total_delta_r2_comparisons = nrow(current_delta),
      stringsAsFactors = FALSE
    )
  }

  main_table <- do.call(rbind, main_rows)
  trait_detail <- do.call(rbind, trait_detail_rows)
  pairwise_overlap <- do.call(rbind, pairwise_rows)
  delta_r2_detail <- do.call(rbind, delta_rows)

  rownames(main_table) <- NULL
  rownames(trait_detail) <- NULL
  rownames(pairwise_overlap) <- NULL
  rownames(delta_r2_detail) <- NULL

  methods <- data.frame(
    item = c(
      "Observation record",
      "Observed species",
      "Union species",
      "Species shared by all traits",
      "Shared fraction of union",
      "Pairwise shared species",
      "Delta R2",
      "Positive Delta R2",
      "Negative Delta R2"
    ),
    definition = c(
      paste0(
        "One non-missing trait value in one row of trait_data after ",
        "restricting species to tree tips."
      ),
      "A species with at least one non-missing record for a trait.",
      "The number of species observed for at least one trait in the set.",
      paste0(
        "The number of species with at least one observation for every ",
        "trait in the set; records do not need to come from the same row."
      ),
      "Species shared by all traits divided by union species.",
      "The number of observed species shared by each pair of traits.",
      "Multi-trait R2 minus the corresponding single-trait R2.",
      "Delta R2 is greater than delta_tolerance.",
      "Delta R2 is less than negative delta_tolerance."
    ),
    stringsAsFactors = FALSE
  )

  tables <- list(
    `Table 1` = main_table,
    `Trait Detail` = trait_detail,
    `Pairwise Overlap` = pairwise_overlap,
    `Delta R2 Detail` = delta_r2_detail,
    Methods = methods
  )

  saved_files <- character(0L)
  if (!is.null(output_dir)) {
    if (!is.character(output_dir) || length(output_dir) != 1L) {
      stop_with("output_dir must be NULL or one directory path.")
    }

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output_dir)) {
      stop_with("Could not create output directory: %s", output_dir)
    }

    csv_paths <- c(
      main_table = file.path(
        output_dir,
        paste0(file_prefix, "_table1.csv")
      ),
      trait_detail = file.path(
        output_dir,
        paste0(file_prefix, "_trait_detail.csv")
      ),
      pairwise_overlap = file.path(
        output_dir,
        paste0(file_prefix, "_pairwise_overlap.csv")
      ),
      delta_r2_detail = file.path(
        output_dir,
        paste0(file_prefix, "_delta_r2_detail.csv")
      )
    )

    write_csv_file(main_table, csv_paths[["main_table"]])
    write_csv_file(trait_detail, csv_paths[["trait_detail"]])
    write_csv_file(pairwise_overlap, csv_paths[["pairwise_overlap"]])
    write_csv_file(delta_r2_detail, csv_paths[["delta_r2_detail"]])
    saved_files <- c(saved_files, csv_paths)

    result_path <- file.path(
      output_dir,
      paste0(file_prefix, ".rds")
    )

    result_object <- list(
      main_table = main_table,
      trait_detail = trait_detail,
      pairwise_overlap = pairwise_overlap,
      delta_r2_detail = delta_r2_detail,
      methods = methods
    )
    saveRDS(result_object, result_path)
    saved_files <- c(saved_files, result_rds = result_path)

    if (isTRUE(save_xlsx)) {
      workbook_path <- file.path(
        output_dir,
        paste0(file_prefix, ".xlsx")
      )

      if (write_excel_file(tables, workbook_path)) {
        saved_files <- c(saved_files, workbook = workbook_path)
      }
    }
  }

  result <- list(
    main_table = main_table,
    trait_detail = trait_detail,
    pairwise_overlap = pairwise_overlap,
    delta_r2_detail = delta_r2_detail,
    methods = methods,
    saved_files = saved_files
  )

  class(result) <- c("results2_summary_tables", class(result))
  result
}

# ------------------------------------------------------------
# Function: Extract selected traits with optional log transformation
# ------------------------------------------------------------

extract_traits <- function(trait_data,
                           trait_names,
                           log_transform = FALSE,
                           scale_transform = FALSE) {

  # Check whether species column exists
  if (!"species" %in% names(trait_data)) {
    stop("The dataset must contain a column named 'species'.")
  }

  # Check whether all selected traits exist
  missing_traits <- setdiff(trait_names, names(trait_data))

  if (length(missing_traits) > 0) {
    stop(paste(
        "The following trait(s) are not found in trait_data:",
        paste(missing_traits, collapse = ", "))
    )
  }

  # Extract species and selected traits
  selected_trait_data <- trait_data %>%
    dplyr::select(species, dplyr::all_of(trait_names))

  # Check whether selected traits are numeric
  non_numeric_traits <- trait_names[
    !vapply(selected_trait_data[trait_names],
            is.numeric,
            logical(1))
  ]

  if (length(non_numeric_traits) > 0) {
    stop(paste(
        "The following trait(s) are not numeric:",
        paste(non_numeric_traits, collapse = ", "))
    )
  }

  # Apply log transformation if requested
  if (log_transform) {

    non_positive_traits <- trait_names[
      vapply(selected_trait_data[trait_names],
             function(x) any(x <= 0, na.rm = TRUE),
             logical(1))
    ]

    if (length(non_positive_traits) > 0) {
      stop(paste(paste0(
            "Log transformation requires positive values. ",
            "The following trait(s) contain zero or negative values:"),
          paste(non_positive_traits, collapse = ", "))
      )
    }

    selected_trait_data <- selected_trait_data %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(trait_names),
          log
        )
      )
  }

  # Apply scale transformation if requested
  if (scale_transform) {

    zero_variance_traits <- trait_names[
      vapply(selected_trait_data[trait_names],
             function(x) {
               x <- x[!is.na(x)]
               length(x) < 2 || stats::sd(x) == 0
             },
             logical(1))
    ]

    if (length(zero_variance_traits) > 0) {
      stop(paste(
        paste0(
          "The following trait(s) cannot be scaled because ",
          "they have fewer than two observed values or zero variance:"
        ),
        paste(zero_variance_traits, collapse = ", "))
      )
    }

    selected_trait_data <- selected_trait_data %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(trait_names),
          ~ as.numeric(scale(.x))
        )
      )
  }

  # Store metadata
  attr(selected_trait_data, "trait_names") <- trait_names
  attr(selected_trait_data, "log_transform") <- log_transform
  attr(selected_trait_data, "scale_transform") <- scale_transform

  # Add class for custom printing
  class(selected_trait_data) <- c("extracted_traits", class(selected_trait_data))

  return(selected_trait_data)
}

# ------------------------------------------------------------
# Print method for extracted trait data
# ------------------------------------------------------------

print.extracted_traits <- function(x, ...) {

  cat("Selected traits: ",
      paste(attr(x, "trait_names"), collapse = ", "),
      "\n",
      sep = "")

  cat("Log transformed: ",
      attr(x, "log_transform"),
      "\n",
      sep = "")

  cat("Scale transformed: ",
      attr(x, "scale_transform"),
      "\n\n",
      sep = "")

  NextMethod()
}


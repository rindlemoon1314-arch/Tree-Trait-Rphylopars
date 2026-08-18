# ------------------------------------------------------------
# Function: Create k-fold CV datasets for trait imputation
# ------------------------------------------------------------

create_kfold_cv <- function(trait_data,
                            k = 5,
                            seed = 123,
                            scale_transform = FALSE) {

  # Check input data
  if (!is.data.frame(trait_data)) {
    stop("trait_data must be a data frame.")
  }

  if (!"species" %in% names(trait_data)) {
    stop("The input dataset must contain a column named 'species'.")
  }

  if (length(k) != 1 ||
      !is.numeric(k) ||
      k < 2 ||
      k %% 1 != 0) {
    stop("k must be a single integer greater than or equal to 2.")
  }

  # Identify trait columns
  trait_names <- setdiff(names(trait_data), "species")

  if (length(trait_names) == 0) {
    stop("No trait columns were found.")
  }

  # Check whether traits are numeric
  non_numeric_traits <- trait_names[
    !vapply(trait_data[trait_names],
            is.numeric,
            logical(1))
  ]

  if (length(non_numeric_traits) > 0) {
    stop(
      "The following trait(s) are not numeric: ",
      paste(non_numeric_traits, collapse = ", ")
    )
  }

  set.seed(seed)

  n <- nrow(trait_data)
  p <- length(trait_names)

  # Store fold membership for every observed trait value
  fold_id <- matrix(NA_integer_,
                    nrow = n,
                    ncol = p,
                    dimnames = list(NULL, trait_names))

  # Independently assign the observed values of each trait to one of the k folds
  for (j in seq_along(trait_names)) {

    trait <- trait_names[j]
    observed_rows <- which(!is.na(trait_data[[trait]]))

    if (length(observed_rows) == 0) {
      warning(sprintf("Trait '%s' contains no observed values.", trait))
      next
    }

    if (length(observed_rows) < k) {
      warning(sprintf(paste0(
            "Trait '%s' has fewer observed values than k; ",
            "some folds will contain no held-out value."),
          trait
        )
      )
    }

    fold_id[observed_rows, j] <- sample(rep(seq_len(k), length.out = length(observed_rows)))
  }

  # Create the CV datasets
  cv_datasets <- vector("list", k)

  for (fold in seq_len(k)) {

    fold_data <- trait_data

    scale_parameters <- vector("list",
                               length(trait_names))
    names(scale_parameters) <- trait_names

    assessment <- vector("list",
                         length(trait_names))

    for (j in seq_along(trait_names)) {

      trait <- trait_names[j]

      # Values held out in this fold
      rows_to_hide <- which(fold_id[, j] == fold)

      # Training values exclude held-out observations
      training_values <- trait_data[[trait]]
      training_values[rows_to_hide] <- NA_real_

      if (scale_transform) {

        # Calculate parameters from training values only
        center <- mean(training_values, na.rm = TRUE)
        scale_value <- stats::sd(training_values,
                                 na.rm = TRUE)

        if (!is.finite(center) ||
            !is.finite(scale_value) ||
            scale_value == 0) {
          stop(sprintf(paste0(
                "Trait '%s' cannot be scaled in fold %d ",
                "because the training values have zero ",
                "variance or insufficient observations."),
              trait,
              fold
            )
          )
        }

        # Apply training parameters to the complete column
        fold_data[[trait]] <- (trait_data[[trait]] - center) / scale_value

        transformed_truth <- (trait_data[[trait]][rows_to_hide] - center) / scale_value

        scale_parameters[[trait]] <- c(center = center,
                                       scale = scale_value)

      } else {

        transformed_truth <- trait_data[[trait]][rows_to_hide]

        scale_parameters[[trait]] <- c(center = 0,
                                       scale = 1)
      }

      # Hide validation values after transformation
      fold_data[[trait]][rows_to_hide] <- NA_real_

      # Save the true held-out values
      assessment[[j]] <- data.frame(row_id = rows_to_hide,
                                    species = trait_data$species[rows_to_hide],
                                    trait = trait,
                                    raw_truth = trait_data[[trait]][rows_to_hide],
                                    transformed_truth = transformed_truth,
                                    stringsAsFactors = FALSE)
    }

    cv_datasets[[fold]] <- list(data = fold_data,
                                assessment = do.call(rbind, assessment),
                                scale_parameters = scale_parameters)
  }

  names(cv_datasets) <- paste0("fold_",
                               seq_len(k))

  # -------------------------------------------------------
  # Create summary of held-out observations
  # -------------------------------------------------------

  fold_summary <- tibble(trait = trait_names,
                         total_observations = as.integer(colSums(!is.na(trait_data[trait_names])))
  )

  # Number of values changed to NA in each fold
  for (fold in seq_len(k)) {

    fold_summary[[paste0("fold_", fold)]] <-
      as.integer(colSums(fold_id == fold,na.rm = TRUE))
  }

  attr(cv_datasets, "fold_id") <- fold_id
  attr(cv_datasets, "scale_transform") <- scale_transform
  attr(cv_datasets, "fold_summary") <- fold_summary

  #  Add class for custom printing
  class(cv_datasets) <- c("kfold_cv_datasets",
                          class(cv_datasets))

  return(cv_datasets)
}


# ---------------------------------------------------------
# Print method for k-fold CV datasets
# ---------------------------------------------------------

print.kfold_cv_datasets <- function(x, ...) {

  cat("K-fold CV datasets\n")
  cat("Number of folds:",
      length(x),
      "\n")
  print(attr(x, "fold_summary"))
  invisible(x)
}


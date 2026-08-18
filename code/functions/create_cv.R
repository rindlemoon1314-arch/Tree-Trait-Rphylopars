# ------------------------------------------------------------
# Function: create cross-validation datasets for one or multiple traits
# ------------------------------------------------------------

create_cross_validation <- function(trait_data, test_proportion = 0.25, seed = 123) {

  # Check whether the dataset contains a species column
  if (!"species" %in% names(trait_data)) {
    stop("The input dataset must contain a column named 'species'.")
  }

  # Identify trait columns
  trait_names <- names(trait_data)[names(trait_data) != "species"]

  # Create the cross-validation dataset
  trait_cv <- trait_data

  # Create an empty list to store test sets
  testset_list <- list()

  # Set random seed for reproducibility
  set.seed(seed)

  # Loop through each trait
  for (trait in trait_names) {

    # Find rows with observed values for this trait
    observed_rows <- which(!is.na(trait_data[[trait]]))

    # Randomly select test rows
    test_rows <- sample(
      observed_rows,
      size = round(length(observed_rows) * test_proportion),
      replace = FALSE
    )

    # Save the test set in long format
    testset_trait <- trait_data[test_rows, ] %>%
      dplyr::select(species, dplyr::all_of(trait)) %>%
      dplyr::rename(observed = dplyr::all_of(trait)) %>%
      dplyr::mutate(trait = trait) %>%
      dplyr::select(species, trait, observed)

    # Set the selected test values to NA in the CV dataset
    trait_cv[[trait]][test_rows] <- NA

    # Store the test set
    testset_list[[trait]] <- testset_trait
  }

  # Combine all test sets
  testset <- dplyr::bind_rows(testset_list)

  # Return the CV dataset and test set
  return(
    list(
      trait_cv = trait_cv,
      testset = testset
    )
  )
}

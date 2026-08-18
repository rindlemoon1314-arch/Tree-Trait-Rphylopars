# ------------------------------------------------------------
# Function: Plot raw and log-transformed trait distributions
# ------------------------------------------------------------

plot_trait_histograms <- function(trait_data,
                                  trait_names,
                                  bins = 40) {

  # Check input data
  if (!is.data.frame(trait_data)) {
    stop("trait_data must be a data frame.")
  }

  # Check whether selected traits exist
  missing_traits <- setdiff(trait_names,
                            names(trait_data))

  if (length(missing_traits) > 0) {
    stop(paste("The following trait(s) were not found:",
         paste(missing_traits,
               collapse = ", ")))
  }

  # Check whether selected traits are numeric
  non_numeric_traits <- trait_names[!vapply(
    trait_data[trait_names],
    is.numeric,
    logical(1))]

  if (length(non_numeric_traits) > 0) {
    stop(paste("The following trait(s) are not numeric:",
         paste(non_numeric_traits,
               collapse = ", "))
    )
  }

  # Check whether log transformation is possible
  non_positive_traits <- trait_names[vapply(
    trait_data[trait_names],
    function(x) {
        any(x <= 0,
            na.rm = TRUE)
    },
    logical(1)
    )
  ]

  if (length(non_positive_traits) > 0) {
    stop(paste(paste0(
      "Log transformation requires positive values. ",
      "The following trait(s) contain zero or negative values:"),
      paste(non_positive_traits,
            collapse = ", ")
      )
    )
  }

  # Create an empty list for plots
  plot_list <- vector("list", length(trait_names))

  names(plot_list) <- trait_names

  # Create two plots for each trait
  for (trait in trait_names) {

    values <- trait_data[[trait]]

    # Remove missing values
    values <- values[!is.na(values)]

    # Prepare raw and log-transformed data
    raw_data <- data.frame(value = values)
    log_data <- data.frame(value = log(values))

    # --------------------------------------------------------
    # Raw-value histogram with density curve
    # --------------------------------------------------------

    raw_plot <- ggplot(raw_data, aes(x = value)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = bins,
                     fill = "#3B82B4",
                     colour = "white",
                     alpha = 0.7) +
      geom_density(colour = "#D73027",
                   linewidth = 1) +
      labs(x = trait,
           y = "Density",
           title = paste0(trait, ": raw values")) +
      theme_bw()

    # --------------------------------------------------------
    # Log-transformed histogram with density curve
    # --------------------------------------------------------

    log_plot <- ggplot(log_data, aes(x = value)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = bins,
                     fill = "#4DAF4A",
                     colour = "white",
                     alpha = 0.7) +
      geom_density(colour = "#D73027",
                   linewidth = 1) +
      labs(x = paste0("log(", trait, ")"),
           y = "Density",
           title = paste0(trait, ": log-transformed values")) +
      theme_bw()

    # Store both plots
    plot_list[[trait]] <- list(raw = raw_plot, log = log_plot)
  }

  return(plot_list)
}

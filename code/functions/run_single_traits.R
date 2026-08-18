# ------------------------------------------------------------
# Run k-fold phylopars analyses for one or more individual traits
#
# Sourcing this file only defines functions. It does not start an
# analysis automatically.
# ------------------------------------------------------------


# Default analysis configurations used in the project.
single_trait_default_configurations <- function() {

  data.frame(
    config_id = c(
      "BM_FTF",
      "BM_TTF",
      "lambda_FTF",
      "lambda_TTF"
    ),
    model = c(
      "BM",
      "BM",
      "lambda",
      "lambda"
    ),
    pheno_error = c(
      FALSE,
      TRUE,
      FALSE,
      TRUE
    ),
    phylo_correlated = TRUE,
    pheno_correlated = FALSE,
    stringsAsFactors = FALSE
  )
}


# Keep the original fold assignments while retaining only one trait.
make_single_trait_cv <- function(cv_datasets, trait) {

  output <- lapply(
    cv_datasets,
    function(fold) {

      list(
        data = fold$data[
          ,
          c("species", trait),
          drop = FALSE
        ],
        assessment = fold$assessment[
          fold$assessment$trait == trait,
          ,
          drop = FALSE
        ],
        scale_parameters = fold$scale_parameters[trait]
      )
    }
  )

  names(output) <- names(cv_datasets)

  original_fold_id <- attr(cv_datasets, "fold_id")

  if (!is.null(original_fold_id)) {
    attr(output, "fold_id") <- original_fold_id[
      ,
      trait,
      drop = FALSE
    ]
  }

  attr(output, "scale_transform") <- attr(
    cv_datasets,
    "scale_transform"
  )

  original_summary <- attr(cv_datasets, "fold_summary")

  if (!is.null(original_summary)) {
    attr(output, "fold_summary") <- original_summary[
      original_summary$trait == trait,
      ,
      drop = FALSE
    ]
  }

  class(output) <- class(cv_datasets)
  output
}


# Save through a temporary file so that an interrupted write is less
# likely to leave a corrupt final result file.
save_rds_atomic <- function(object, path, compress = TRUE) {

  temporary_path <- tempfile(
    pattern = paste0(".", basename(path), "_"),
    tmpdir = dirname(path)
  )

  on.exit(
    unlink(temporary_path),
    add = TRUE
  )

  saveRDS(
    object,
    temporary_path,
    compress = compress
  )

  if (file.exists(path) && !file.remove(path)) {
    stop("Could not replace the existing result file: ", path)
  }

  if (!file.rename(temporary_path, path)) {
    stop("Could not move the completed result into place: ", path)
  }

  invisible(path)
}


# ------------------------------------------------------------
# Main function
# ------------------------------------------------------------

run_single_traits <- function(
    cv_datasets,
    tree,
    traits,
    result_dir,
    checkpoint_dir = file.path(result_dir, "checkpoints"),
    configurations = single_trait_default_configurations(),
    runner = run_kfold_phylopars_parallel,
    n_cores = 2L,
    resume = TRUE,
    retry_failed = FALSE,
    worker_stop_grace = 0.5,
    result_components = c("anc_recon", "anc_var"),
    overwrite = FALSE,
    file_prefix = "PHY_kfoldcv",
    compress = TRUE
) {

  # ----------------------------------------------------------
  # Validate inputs
  # ----------------------------------------------------------

  if (!is.list(cv_datasets) || length(cv_datasets) == 0L) {
    stop("cv_datasets must be a non-empty list of CV folds.")
  }

  valid_folds <- vapply(
    cv_datasets,
    function(fold) {
      is.list(fold) &&
        is.data.frame(fold$data) &&
        is.data.frame(fold$assessment) &&
        is.list(fold$scale_parameters) &&
        "species" %in% names(fold$data) &&
        "trait" %in% names(fold$assessment)
    },
    logical(1)
  )

  if (!all(valid_folds)) {
    stop(
      "Every CV fold must contain valid data, assessment, and ",
      "scale_parameters components."
    )
  }

  if (
    !is.character(traits) ||
    length(traits) == 0L ||
    anyNA(traits) ||
    any(!nzchar(traits))
  ) {
    stop("traits must be a non-empty character vector.")
  }

  traits <- unique(traits)

  missing_trait_folds <- vapply(
    traits,
    function(trait) {
      any(
        !vapply(
          cv_datasets,
          function(fold) {
            trait %in% names(fold$data) &&
              trait %in% fold$assessment$trait &&
              trait %in% names(fold$scale_parameters)
          },
          logical(1)
        )
      )
    },
    logical(1)
  )

  if (any(missing_trait_folds)) {
    stop(
      "These traits are incomplete or absent in cv_datasets: ",
      paste(traits[missing_trait_folds], collapse = ", ")
    )
  }

  if (
    length(result_dir) != 1L ||
    is.na(result_dir) ||
    !nzchar(result_dir)
  ) {
    stop("result_dir must be one non-empty path.")
  }

  if (
    length(checkpoint_dir) != 1L ||
    is.na(checkpoint_dir) ||
    !nzchar(checkpoint_dir)
  ) {
    stop("checkpoint_dir must be one non-empty path.")
  }

  if (!is.function(runner)) {
    stop("runner must be a function.")
  }

  required_runner_arguments <- c(
    "cv_datasets",
    "tree",
    "model",
    "pheno_error",
    "phylo_correlated",
    "pheno_correlated",
    "usezscores",
    "EM_verbose",
    "optim_verbose",
    "n_cores",
    "checkpoint_dir",
    "resume",
    "retry_failed",
    "worker_stop_grace",
    "result_components"
  )

  runner_arguments <- names(formals(runner))

  if (!"..." %in% runner_arguments) {
    missing_runner_arguments <- setdiff(
      required_runner_arguments,
      runner_arguments
    )

    if (length(missing_runner_arguments) > 0L) {
      stop(
        "runner is missing these required arguments: ",
        paste(missing_runner_arguments, collapse = ", ")
      )
    }
  }

  required_configuration_columns <- c(
    "config_id",
    "model",
    "pheno_error",
    "phylo_correlated",
    "pheno_correlated"
  )

  if (!is.data.frame(configurations)) {
    stop("configurations must be a data frame.")
  }

  missing_configuration_columns <- setdiff(
    required_configuration_columns,
    names(configurations)
  )

  if (length(missing_configuration_columns) > 0L) {
    stop(
      "configurations is missing these columns: ",
      paste(missing_configuration_columns, collapse = ", ")
    )
  }

  if (nrow(configurations) == 0L) {
    stop("configurations must contain at least one row.")
  }

  if (
    anyNA(configurations$config_id) ||
    any(!nzchar(configurations$config_id)) ||
    anyDuplicated(configurations$config_id)
  ) {
    stop("configurations$config_id must contain unique names.")
  }

  logical_configuration_columns <- c(
    "pheno_error",
    "phylo_correlated",
    "pheno_correlated"
  )

  invalid_logical_columns <- logical_configuration_columns[
    !vapply(
      configurations[logical_configuration_columns],
      function(x) is.logical(x) && !anyNA(x),
      logical(1)
    )
  ]

  if (length(invalid_logical_columns) > 0L) {
    stop(
      "These configuration columns must contain logical values: ",
      paste(invalid_logical_columns, collapse = ", ")
    )
  }

  if (
    length(n_cores) != 1L ||
    !is.numeric(n_cores) ||
    !is.finite(n_cores) ||
    n_cores < 1
  ) {
    stop("n_cores must be a positive integer.")
  }

  n_cores <- as.integer(n_cores)

  logical_options <- list(
    resume = resume,
    retry_failed = retry_failed,
    overwrite = overwrite
  )

  invalid_logical_options <- names(logical_options)[
    !vapply(
      logical_options,
      function(x) {
        length(x) == 1L && is.logical(x) && !is.na(x)
      },
      logical(1)
    )
  ]

  if (length(invalid_logical_options) > 0L) {
    stop(
      "These arguments must be TRUE or FALSE: ",
      paste(invalid_logical_options, collapse = ", ")
    )
  }

  safe_name <- function(x) {
    gsub("[^A-Za-z0-9_]+", "_", x)
  }

  trait_ids <- safe_name(traits)
  configuration_ids <- safe_name(configurations$config_id)

  if (anyDuplicated(trait_ids)) {
    stop("Some trait names become identical after filename cleaning.")
  }

  if (anyDuplicated(configuration_ids)) {
    stop("Some config_id values become identical after filename cleaning.")
  }

  dir.create(
    result_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    checkpoint_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!dir.exists(result_dir)) {
    stop("Could not create result_dir: ", result_dir)
  }

  if (!dir.exists(checkpoint_dir)) {
    stop("Could not create checkpoint_dir: ", checkpoint_dir)
  }

  # ----------------------------------------------------------
  # Build and execute the run plan
  # ----------------------------------------------------------

  run_plan <- expand.grid(
    trait_index = seq_along(traits),
    configuration_index = seq_len(nrow(configurations)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  run_plan <- run_plan[
    order(
      run_plan$trait_index,
      run_plan$configuration_index
    ),
    ,
    drop = FALSE
  ]

  run_summary <- vector(
    "list",
    nrow(run_plan)
  )

  for (plan_index in seq_len(nrow(run_plan))) {

    trait_index <- run_plan$trait_index[[plan_index]]
    configuration_index <- run_plan$configuration_index[[plan_index]]

    trait <- traits[[trait_index]]
    trait_id <- trait_ids[[trait_index]]

    configuration <- configurations[
      configuration_index,
      ,
      drop = FALSE
    ]

    configuration_id <- configuration_ids[[configuration_index]]
    run_id <- paste(trait_id, configuration_id, sep = "__")

    current_checkpoint_dir <- file.path(
      checkpoint_dir,
      paste0("kfold_", trait_id, "_", configuration_id)
    )

    result_path <- file.path(
      result_dir,
      paste0(
        safe_name(file_prefix),
        "_",
        trait_id,
        "_",
        configuration_id,
        ".rds"
      )
    )

    message("\n============================================================")
    message(
      "Run ",
      plan_index,
      " of ",
      nrow(run_plan),
      ": ",
      run_id
    )
    message("============================================================")

    if (file.exists(result_path) && !overwrite) {
      message("Already saved; skipping: ", result_path)

      run_summary[[plan_index]] <- data.frame(
        trait = trait,
        configuration = configuration$config_id[[1L]],
        model = configuration$model[[1L]],
        state = "already_saved",
        successful_folds = NA_integer_,
        total_folds = length(cv_datasets),
        result_path = result_path,
        checkpoint_path = current_checkpoint_dir,
        stringsAsFactors = FALSE
      )

      next
    }

    single_trait_cv <- make_single_trait_cv(
      cv_datasets,
      trait
    )

    fit <- runner(
      cv_datasets = single_trait_cv,
      tree = tree,
      model = configuration$model[[1L]],
      pheno_error = configuration$pheno_error[[1L]],
      phylo_correlated = configuration$phylo_correlated[[1L]],
      pheno_correlated = configuration$pheno_correlated[[1L]],
      usezscores = FALSE,
      EM_verbose = FALSE,
      optim_verbose = FALSE,
      n_cores = n_cores,
      checkpoint_dir = current_checkpoint_dir,
      resume = resume,
      retry_failed = retry_failed,
      worker_stop_grace = worker_stop_grace,
      result_components = result_components
    )

    if (isTRUE(attr(fit, "interrupted"))) {
      rm(single_trait_cv, fit)
      gc(verbose = FALSE)

      stop(
        "Run ",
        run_id,
        " was interrupted. Its checkpoints were retained. ",
        "Call run_single_traits() again with the same paths to resume."
      )
    }

    status_table <- attr(fit, "status")

    if (is.null(status_table) || !"status" %in% names(status_table)) {
      stop("runner returned a result without a valid status table.")
    }

    successful_folds <- sum(
      status_table$status == "success",
      na.rm = TRUE
    )

    save_rds_atomic(
      fit,
      result_path,
      compress = compress
    )

    result_state <- if (successful_folds == nrow(status_table)) {
      "saved_complete"
    } else {
      "saved_with_failed_folds"
    }

    run_summary[[plan_index]] <- data.frame(
      trait = trait,
      configuration = configuration$config_id[[1L]],
      model = configuration$model[[1L]],
      state = result_state,
      successful_folds = successful_folds,
      total_folds = nrow(status_table),
      result_path = result_path,
      checkpoint_path = current_checkpoint_dir,
      stringsAsFactors = FALSE
    )

    message(
      "Saved: ",
      result_path,
      " (",
      successful_folds,
      "/",
      nrow(status_table),
      " successful folds)"
    )

    rm(single_trait_cv, fit, status_table)
    gc(verbose = FALSE)
  }

  run_summary <- do.call(
    rbind,
    run_summary
  )

  rownames(run_summary) <- NULL

  message(
    "\nBatch finished. Result files are in:\n",
    normalizePath(
      result_dir,
      winslash = "/",
      mustWork = FALSE
    )
  )

  invisible(run_summary)
}


# ------------------------------------------------------------
# Function: Run phylopars for k-fold CV datasets in parallel
#           with per-fold checkpoints and interrupt recovery
# ------------------------------------------------------------

run_kfold_phylopars_parallel <- function(
    cv_datasets,
    tree,
    model = "BM",
    pheno_error = TRUE,
    phylo_correlated = TRUE,
    phylocov_start = NULL,
    pheno_correlated = TRUE,
    usezscores = TRUE,
    EM_verbose = TRUE,
    optim_verbose = TRUE,
    n_cores = max(1L, parallel::detectCores() - 2L),
    checkpoint_dir = NULL,
    resume = TRUE,
    retry_failed = FALSE,
    worker_stop_grace = 0.5,
    result_components = c("anc_recon", "anc_var")
) {

  # ----------------------------------------------------------
  # Input validation
  # ----------------------------------------------------------

  if (!is.list(cv_datasets) || length(cv_datasets) == 0L) {
    stop("cv_datasets must be a non-empty list.")
  }

  valid_folds <- vapply(
    cv_datasets,
    function(x) {
      is.list(x) &&
        "data" %in% names(x) &&
        is.data.frame(x$data)
    },
    logical(1)
  )

  if (!all(valid_folds)) {
    stop(
      "Each fold must be a list containing a data frame named 'data'."
    )
  }

  if (
    length(usezscores) != 1L ||
    !is.logical(usezscores) ||
    is.na(usezscores)
  ) {
    stop("usezscores must be TRUE or FALSE.")
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

  if (
    length(resume) != 1L ||
    !is.logical(resume) ||
    is.na(resume)
  ) {
    stop("resume must be TRUE or FALSE.")
  }

  if (
    length(retry_failed) != 1L ||
    !is.logical(retry_failed) ||
    is.na(retry_failed)
  ) {
    stop("retry_failed must be TRUE or FALSE.")
  }

  if (
    length(worker_stop_grace) != 1L ||
    !is.numeric(worker_stop_grace) ||
    !is.finite(worker_stop_grace) ||
    worker_stop_grace < 0
  ) {
    stop("worker_stop_grace must be a non-negative number of seconds.")
  }

  if (
    !is.null(result_components) &&
    (
      !is.character(result_components) ||
      length(result_components) == 0L ||
      anyNA(result_components) ||
      any(!nzchar(result_components))
    )
  ) {
    stop(
      paste0(
        "result_components must be NULL or a non-empty character vector."
      )
    )
  }

  n_folds <- length(cv_datasets)

  fold_names <- names(cv_datasets)

  if (is.null(fold_names)) {
    fold_names <- paste0("fold_", seq_len(n_folds))
  } else {
    missing_names <- is.na(fold_names) | fold_names == ""
    fold_names[missing_names] <- paste0("fold_", which(missing_names))
  }

  # ----------------------------------------------------------
  # Starting covariance matrices
  # ----------------------------------------------------------

  if (is.null(phylocov_start)) {
    fold_start_cov <- rep(list(NULL), n_folds)
  } else if (is.matrix(phylocov_start)) {
    fold_start_cov <- rep(list(phylocov_start), n_folds)
  } else if (
    is.list(phylocov_start) &&
    length(phylocov_start) == n_folds
  ) {
    fold_start_cov <- phylocov_start
  } else {
    stop(
      paste0(
        "phylocov_start must be NULL, a matrix, ",
        "or a list with one element per fold."
      )
    )
  }

  names(fold_start_cov) <- fold_names

  # ----------------------------------------------------------
  # Checkpoint directory and run manifest
  # ----------------------------------------------------------

  temporary_checkpoint_dir <- is.null(checkpoint_dir)

  if (temporary_checkpoint_dir) {
    checkpoint_dir <- tempfile(
      pattern = "kfold_phylopars_checkpoints_"
    )
  }

  checkpoint_dir <- normalizePath(
    checkpoint_dir,
    winslash = "/",
    mustWork = FALSE
  )

  if (!dir.exists(checkpoint_dir)) {
    dir.create(
      checkpoint_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }

  if (!dir.exists(checkpoint_dir)) {
    stop("Could not create checkpoint directory: ", checkpoint_dir)
  }

  message("Checkpoint directory: ", checkpoint_dir)

  if (temporary_checkpoint_dir) {
    message(
      paste0(
        "For recovery after restarting R, pass an explicit persistent ",
        "checkpoint_dir."
      )
    )
  }

  checkpoint_width <- max(4L, nchar(as.character(n_folds)))

  checkpoint_files <- file.path(
    checkpoint_dir,
    sprintf(
      paste0("fold_%0", checkpoint_width, "d.rds"),
      seq_len(n_folds)
    )
  )

  manifest_file <- file.path(checkpoint_dir, "manifest.rds")

  object_md5 <- function(x) {
    signature_file <- tempfile(
      pattern = "kfold_signature_",
      fileext = ".rds"
    )

    on.exit(unlink(signature_file), add = TRUE)

    saveRDS(
      x,
      signature_file,
      version = 3,
      compress = FALSE
    )

    unname(tools::md5sum(signature_file))
  }

  signature_object <- list(
    cv_datasets = cv_datasets,
    tree = tree,
    model = model,
    pheno_error = pheno_error,
    phylo_correlated = phylo_correlated,
    phylocov_start = fold_start_cov,
    pheno_correlated = pheno_correlated,
    usezscores = usezscores
  )

  message("Validating run signature...")

  run_signature <- object_md5(signature_object)

  if (file.exists(manifest_file)) {

    if (!resume) {
      stop(
        paste0(
          "The checkpoint directory already contains a run. ",
          "Set resume = TRUE or use a new checkpoint_dir."
        )
      )
    }

    manifest <- tryCatch(
      readRDS(manifest_file),
      error = function(e) NULL
    )

    if (
      is.null(manifest) ||
      !identical(manifest$signature, run_signature)
    ) {
      stop(
        paste0(
          "The existing checkpoint directory belongs to a different ",
          "dataset, tree, or model configuration. Use a different ",
          "checkpoint_dir."
        )
      )
    }

  } else {

    existing_files <- list.files(
      checkpoint_dir,
      all.files = TRUE,
      no.. = TRUE
    )

    if (length(existing_files) > 0L) {
      stop(
        paste0(
          "checkpoint_dir is not empty and does not contain a valid ",
          "manifest. Use an empty directory."
        )
      )
    }

    manifest <- list(
      version = 2L,
      signature = run_signature,
      fold_names = fold_names,
      n_folds = n_folds,
      model = model,
      created_at = Sys.time()
    )

    saveRDS(manifest, manifest_file)
  }

  # ----------------------------------------------------------
  # Checkpoint helpers
  # ----------------------------------------------------------

  read_fold_checkpoint <- function(i) {

    path <- checkpoint_files[[i]]

    if (!file.exists(path)) {
      return(NULL)
    }

    checkpoint <- tryCatch(
      readRDS(path),
      error = function(e) NULL,
      interrupt = function(e) NULL
    )

    if (
      is.null(checkpoint) ||
      !is.list(checkpoint) ||
      !identical(checkpoint$index, i) ||
      length(checkpoint$status) != 1L ||
      is.na(checkpoint$status) ||
      !checkpoint$status %in% c("success", "failed")
    ) {
      return(NULL)
    }

    checkpoint
  }

  checkpoint_exists <- file.exists(checkpoint_files)

  if (!retry_failed) {
    # Checkpoints are written to a temporary file and atomically renamed.
    # Therefore, an existing final .rds file represents a completed fold.
    # Avoid loading potentially very large fit objects merely to inspect
    # their status.
    checkpoint_is_complete <- checkpoint_exists
  } else {
    message(
      paste0(
        "retry_failed = TRUE: inspecting existing checkpoint statuses..."
      )
    )

    checkpoint_is_complete <- vapply(
      seq_len(n_folds),
      function(i) {
        if (!checkpoint_exists[[i]]) {
          return(FALSE)
        }

        checkpoint <- read_fold_checkpoint(i)
        is_success <- !is.null(checkpoint) &&
          identical(checkpoint$status, "success")

        rm(checkpoint)
        gc(verbose = FALSE)

        is_success
      },
      logical(1)
    )
  }

  pending_folds <- which(!checkpoint_is_complete)

  # Failed or corrupt checkpoint files are removed only when that fold
  # has explicitly been selected for another attempt.
  retry_files <- checkpoint_files[pending_folds]
  retry_files <- retry_files[file.exists(retry_files)]

  if (length(retry_files) > 0L) {
    unlink(retry_files)
  }

  # ----------------------------------------------------------
  # Worker cleanup helpers
  # ----------------------------------------------------------

  cl <- NULL
  pb <- NULL
  worker_pids <- integer(0)
  workers_terminated <- FALSE
  interrupted <- FALSE
  run_problem <- NULL

  pid_is_alive <- function(pid) {
    isTRUE(
      tryCatch(
        tools::pskill(pid, signal = 0L),
        error = function(e) FALSE,
        interrupt = function(e) FALSE
      )
    )
  }

  terminate_worker_processes <- function(pids) {

    pids <- unique(as.integer(pids))
    pids <- pids[is.finite(pids) & pids > 0L]

    if (length(pids) == 0L) {
      return(invisible(NULL))
    }

    # Send SIGTERM only to PIDs obtained directly from this cluster.
    for (pid in pids) {
      tryCatch(
        tools::pskill(pid, signal = tools::SIGTERM),
        error = function(e) NULL,
        interrupt = function(e) NULL
      )
    }

    deadline <- Sys.time() + worker_stop_grace

    repeat {
      alive <- vapply(pids, pid_is_alive, logical(1))

      if (!any(alive) || Sys.time() >= deadline) {
        break
      }

      tryCatch(
        Sys.sleep(0.05),
        interrupt = function(e) NULL
      )
    }

    # A worker inside compiled optimization code may ignore/delay SIGTERM.
    # SIGKILL is restricted to the exact worker PIDs recorded above.
    alive <- vapply(pids, pid_is_alive, logical(1))

    if (any(alive)) {
      for (pid in pids[alive]) {
        tryCatch(
          tools::pskill(pid, signal = tools::SIGKILL),
          error = function(e) NULL,
          interrupt = function(e) NULL
        )
      }
    }

    invisible(NULL)
  }

  abort_workers <- function() {

    if (!workers_terminated) {
      terminate_worker_processes(worker_pids)
      workers_terminated <<- TRUE
    }

    invisible(NULL)
  }

  safe_stop_cluster <- function(cluster) {

    if (is.null(cluster)) {
      return(invisible(NULL))
    }

    tryCatch(
      parallel::stopCluster(cluster),
      error = function(e) NULL,
      interrupt = function(e) NULL
    )

    invisible(NULL)
  }

  safe_close_progress <- function(progress_bar) {

    if (is.null(progress_bar)) {
      return(invisible(NULL))
    }

    tryCatch(
      close(progress_bar),
      error = function(e) NULL,
      interrupt = function(e) NULL
    )

    invisible(NULL)
  }

  # If an unexpected condition escapes the parallel section, terminate
  # only this function's workers before leaving the function.
  on.exit(
    {
      if (!is.null(cl)) {
        abort_workers()
        safe_stop_cluster(cl)
        cl <- NULL
      }

      if (!is.null(pb)) {
        safe_close_progress(pb)
        pb <- NULL
      }
    },
    add = TRUE
  )

  # ----------------------------------------------------------
  # Run pending folds
  # ----------------------------------------------------------

  if (length(pending_folds) > 0L) {

    n_workers <- min(n_cores, length(pending_folds))

    cl <- parallel::makeCluster(n_workers)

    # Record exact worker PIDs before starting any expensive work.
    worker_pids <- unlist(
      parallel::clusterCall(cl, Sys.getpid),
      use.names = FALSE
    )

    doSNOW::registerDoSNOW(cl)

    completed_before_run <- sum(checkpoint_is_complete)

    pb <- utils::txtProgressBar(
      min = 0,
      max = n_folds,
      initial = completed_before_run,
      style = 3
    )

    progress_function <- function(n) {
      utils::setTxtProgressBar(
        pb,
        completed_before_run + n
      )
    }

    opts <- list(progress = progress_function)

    worker_summary <- tryCatch(

      foreach::`%dopar%`(
        foreach::foreach(
          i = pending_folds,
          .packages = "Rphylopars",
          .options.snow = opts,
          .inorder = FALSE
        ),
        {
          fold_data <- cv_datasets[[i]]$data
          start_cov <- fold_start_cov[[i]]

          fold_checkpoint <- tryCatch(
            {
              fit_args <- list(
                trait_data = fold_data,
                tree = tree,
                model = model,
                pheno_error = pheno_error,
                phylo_correlated = phylo_correlated,
                pheno_correlated = pheno_correlated,
                usezscores = usezscores,
                EM_verbose = EM_verbose,
                optim_verbose = optim_verbose
              )

              if (!is.null(start_cov)) {
                fit_args$phylocov_start <- start_cov
              }

              fit <- do.call(
                Rphylopars::phylopars,
                fit_args
              )

              list(
                index = i,
                fold = fold_names[[i]],
                result = fit,
                status = "success",
                error_message = NA_character_,
                completed_at = Sys.time()
              )
            },
            error = function(e) {
              list(
                index = i,
                fold = fold_names[[i]],
                result = NULL,
                status = "failed",
                error_message = conditionMessage(e),
                completed_at = Sys.time()
              )
            }
          )

          checkpoint_path <- checkpoint_files[[i]]

          temporary_path <- paste0(
            checkpoint_path,
            ".",
            Sys.getpid(),
            ".tmp"
          )

          # Write-then-rename prevents a half-written file from being read
          # as a completed fold after interruption.
          saveRDS(
            fold_checkpoint,
            temporary_path,
            version = 3
          )

          saved <- file.rename(
            temporary_path,
            checkpoint_path
          )

          if (!saved) {
            unlink(temporary_path)

            stop(
              "Could not save checkpoint for fold ",
              fold_names[[i]],
              "."
            )
          }

          # Do not transfer the potentially large fit back through the
          # socket. The main process reloads it from the checkpoint.
          list(
            index = i,
            status = fold_checkpoint$status
          )
        }
      ),

      interrupt = function(e) {
        interrupted <<- TRUE

        msg <- conditionMessage(e)

        if (!nzchar(msg)) {
          msg <- "Computation interrupted by the user."
        }

        run_problem <<- msg

        # Esc reaches the master process, not PSOCK workers. Terminate
        # the exact recorded worker PIDs so CPU use stops immediately.
        abort_workers()

        NULL
      },

      error = function(e) {
        run_problem <<- conditionMessage(e)

        # Infrastructure/socket errors must not leave orphan workers.
        abort_workers()

        NULL
      }
    )

    if (!interrupted && is.null(run_problem)) {
      safe_stop_cluster(cl)
    } else {
      abort_workers()
      safe_stop_cluster(cl)
    }

    cl <- NULL
    worker_pids <- integer(0)

    safe_close_progress(pb)
    pb <- NULL

    if (interrupted) {
      message("Run interrupted; returning all completed checkpoints.")
    } else if (!is.null(run_problem)) {
      message(
        "Parallel run stopped; returning all completed checkpoints: ",
        run_problem
      )
    }

  } else {
    message("All folds already have checkpoints; nothing to run.")
  }

  # ----------------------------------------------------------
  # Assemble a full-length result after success or interruption
  # ----------------------------------------------------------

  # Once an interrupt has been caught and the workers have been stopped,
  # protect checkpoint recovery from a queued/second interrupt. Otherwise
  # RStudio can enter Browse[1]> during readRDS(), preventing the assignment
  # on the left-hand side of the function call from being completed.
  if (is.null(result_components)) {
    message(
      paste0(
        "Loading full checkpoint results. This can require substantial RAM..."
      )
    )
  } else {
    message(
      "Loading selected result components: ",
      paste(result_components, collapse = ", "),
      "..."
    )
  }

  final_checkpoints <- suspendInterrupts(
    {
      load_pb <- utils::txtProgressBar(
        min = 0,
        max = n_folds,
        style = 3
      )

      loaded <- lapply(
        seq_len(n_folds),
        function(i) {
          checkpoint <- read_fold_checkpoint(i)

          if (
            !is.null(checkpoint) &&
            !is.null(checkpoint$result) &&
            !is.null(result_components)
          ) {
            missing_components <- setdiff(
              result_components,
              names(checkpoint$result)
            )

            if (length(missing_components) > 0L) {
              checkpoint <- NULL
            } else {
              checkpoint$result <- checkpoint$result[
                result_components
              ]
            }
          }

          utils::setTxtProgressBar(load_pb, i)
          gc(verbose = FALSE)
          checkpoint
        }
      )

      close(load_pb)
      loaded
    }
  )

  missing_message <- if (interrupted) {
    "Fold was not completed because the run was interrupted."
  } else if (!is.null(run_problem)) {
    paste0(
      "Fold was not completed because the parallel run stopped: ",
      run_problem
    )
  } else {
    "Fold did not produce a checkpoint."
  }

  final_checkpoints <- Map(
    function(x, i) {
      if (!is.null(x)) {
        return(x)
      }

      list(
        index = i,
        fold = fold_names[[i]],
        result = NULL,
        status = "failed",
        error_message = missing_message,
        completed_at = as.POSIXct(NA)
      )
    },
    final_checkpoints,
    seq_len(n_folds)
  )

  fold_results <- lapply(
    final_checkpoints,
    `[[`,
    "result"
  )

  names(fold_results) <- fold_names

  status_table <- tibble::tibble(
    fold = fold_names,
    model = rep(as.character(model)[1], n_folds),
    status = vapply(
      final_checkpoints,
      `[[`,
      character(1),
      "status"
    ),
    error_message = vapply(
      final_checkpoints,
      function(x) {
        if (is.null(x$error_message)) {
          NA_character_
        } else {
          as.character(x$error_message)[1]
        }
      },
      character(1)
    ),
    completed_at = as.POSIXct(
      vapply(
        final_checkpoints,
        function(x) as.numeric(x$completed_at),
        numeric(1)
      ),
      origin = "1970-01-01"
    )
  )

  assessments <- lapply(
    cv_datasets,
    function(x) x$assessment
  )

  names(assessments) <- fold_names

  scale_parameters <- lapply(
    cv_datasets,
    function(x) x$scale_parameters
  )

  names(scale_parameters) <- fold_names

  attr(fold_results, "status") <- status_table
  attr(fold_results, "assessment") <- assessments
  attr(fold_results, "scale_parameters") <- scale_parameters
  attr(fold_results, "phylocov_start") <- fold_start_cov
  attr(fold_results, "usezscores") <- usezscores
  attr(fold_results, "model") <- model
  attr(fold_results, "checkpoint_dir") <- checkpoint_dir
  attr(fold_results, "interrupted") <- interrupted
  attr(fold_results, "run_message") <- run_problem
  attr(fold_results, "result_components") <- result_components

  class(fold_results) <- c(
    "kfold_phylopars_results",
    class(fold_results)
  )

  fold_results
}


# ------------------------------------------------------------
# Optional print method
# ------------------------------------------------------------

print.kfold_phylopars_results <- function(x, ...) {

  status_table <- attr(x, "status")

  cat("k-fold phylopars results\n")
  cat("Model:", attr(x, "model"), "\n")
  cat("Interrupted:", isTRUE(attr(x, "interrupted")), "\n")
  cat("Checkpoint directory:", attr(x, "checkpoint_dir"), "\n\n")

  if (!is.null(status_table)) {
    print(status_table)
  }

  invisible(x)
}

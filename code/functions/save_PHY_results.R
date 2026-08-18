# ============================================================
# Safely save PHY result objects one at a time
# ============================================================

save_PHY_results <- function(
    object_names,
    output_dir,
    envir = .GlobalEnv,
    compress = "gzip"
) {

  if (
    !is.character(object_names) ||
    length(object_names) == 0L ||
    anyNA(object_names) ||
    any(!nzchar(object_names))
  ) {
    stop("object_names must be a non-empty character vector.")
  }

  object_names <- unique(object_names)

  if (
    length(output_dir) != 1L ||
    !is.character(output_dir) ||
    is.na(output_dir) ||
    !nzchar(output_dir)
  ) {
    stop("output_dir must be one non-empty path.")
  }

  if (!is.environment(envir)) {
    stop("envir must be an environment.")
  }

  # Check that every requested object exists before saving anything.
  missing_objects <- object_names[
    !vapply(
      object_names,
      exists,
      logical(1),
      envir = envir,
      inherits = FALSE
    )
  ]

  if (length(missing_objects) > 0L) {
    stop(
      "The following objects do not exist; nothing was saved: ",
      paste(missing_objects, collapse = ", ")
    )
  }

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!dir.exists(output_dir)) {
    stop("Could not create the output directory: ", output_dir)
  }

  save_summary <- vector(
    "list",
    length(object_names)
  )

  for (i in seq_along(object_names)) {

    object_name <- object_names[[i]]

    final_path <- file.path(
      output_dir,
      paste0(object_name, ".rds")
    )

    # Retrieve the object from the supplied environment.
    object <- get(
      object_name,
      envir = envir,
      inherits = FALSE
    )

    # Record fold completion status.
    status_table <- attr(object, "status")

    if (
      !is.null(status_table) &&
      is.data.frame(status_table) &&
      "status" %in% names(status_table)
    ) {

      successful_folds <- sum(
        status_table$status == "success",
        na.rm = TRUE
      )

      total_folds <- nrow(status_table)

    } else {

      successful_folds <- sum(
        vapply(
          object,
          Negate(is.null),
          logical(1)
        )
      )

      total_folds <- length(object)
    }

    if (file.exists(final_path)) {

      message(
        "[",
        i,
        "/",
        length(object_names),
        "] File already exists; skipping: ",
        final_path
      )

      save_status <- "already_exists"

    } else {

      # First save to a temporary file in the destination directory.
      temporary_path <- tempfile(
        pattern = paste0(".", object_name, "_"),
        tmpdir = output_dir,
        fileext = ".tmp"
      )

      save_completed <- FALSE

      tryCatch(
        {
          message(
            "[",
            i,
            "/",
            length(object_names),
            "] Saving: ",
            object_name
          )

          saveRDS(
            object,
            file = temporary_path,
            compress = compress
          )

          temporary_size <- file.info(
            temporary_path
          )$size

          if (
            is.na(temporary_size) ||
            temporary_size <= 0
          ) {
            stop("The temporary RDS file is empty.")
          }

          # Rename only after serialization has completed.
          if (!file.rename(temporary_path, final_path)) {
            stop(
              "Could not rename the temporary file to its final name."
            )
          }

          save_completed <- TRUE
          save_status <- "saved"

          message("Saved: ", final_path)
        },
        finally = {
          if (
            !save_completed &&
            file.exists(temporary_path)
          ) {
            unlink(temporary_path)
          }
        }
      )
    }

    final_size <- if (file.exists(final_path)) {
      file.info(final_path)$size / 1024^2
    } else {
      NA_real_
    }

    save_summary[[i]] <- data.frame(
      object = object_name,
      successful_folds = successful_folds,
      total_folds = total_folds,
      status = save_status,
      size_MB = round(final_size, 2),
      file = final_path,
      stringsAsFactors = FALSE
    )

    rm(object, status_table)
    gc(verbose = FALSE)
  }

  save_summary <- do.call(
    rbind,
    save_summary
  )

  rownames(save_summary) <- NULL

  invisible(save_summary)
}

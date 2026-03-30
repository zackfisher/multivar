#' Print dynamics matrix with clear formatting
#'
#' Prints a dynamics matrix (or list of matrices) with clear row/column labels
#' and formatting options. Handles matrices of different sizes gracefully.
#' Can also accept a fit object directly (e.g., from cv.multivar()).
#'
#' @param dynamics A numeric matrix, a list of matrices (for multiple subjects),
#'   or a fit object from cv.multivar() (will automatically extract total effects)
#' @param digits Number of decimal places to display (default: 3)
#' @param zero_char Character to display for zero entries (default: ".")
#' @param highlight_nonzero Logical; if TRUE, emphasize non-zero entries (default: TRUE)
#' @param max_print Maximum number of rows/cols to print before truncating (default: 15)
#' @param time_labels Logical; if TRUE, add subscripts like \eqn{V1_t}, \eqn{V1_{t-1}} (default: TRUE)
#' @param subject_label Optional label for the subject/matrix (e.g., "Subject 1")
#' @param period_labels Optional character vector of period labels for TVP models
#'
#' @return Invisibly returns the input (for chaining)
#' @export
#'
#' @examples
#' # Simple 4x4 matrix
#' phi <- matrix(c(0.3, 0, 0.4, 0, 0.4, 0.5, 0, 0.3, 0, 0.3, 0.2, 0, 0, 0, 0, 0.4), 4, 4)
#' print_dynamics(phi)
#'
#' # List of matrices
#' phi_list <- list(phi, phi * 0.8)
#' print_dynamics(phi_list)
#'
#' # From fit object (automatic extraction)
#' # fit <- cv.multivar(object)
#' # print_dynamics(fit)  # Automatically extracts and prints total effects
#'
#' # TVP fit object with period labels
#' # fit_tvp <- cv.multivar(object_tvp)
#' # print_dynamics(fit_tvp, period_labels = c("Pre", "During", "Post"))
print_dynamics <- function(dynamics,
                          digits = 3,
                          zero_char = ".",
                          highlight_nonzero = TRUE,
                          max_print = 15,
                          time_labels = TRUE,
                          subject_label = NULL,
                          period_labels = NULL) {

  # Detect if input is a fit object (has $mats component)
  if (is.list(dynamics) && !is.null(dynamics$mats) && !is.null(dynamics$mats$total)) {
    fit_obj <- dynamics
    total_effects <- fit_obj$mats$total

    # Detect if TVP (total[[1]] is a list, not a matrix)
    is_tvp <- length(total_effects) > 0 &&
              is.list(total_effects[[1]]) &&
              !is.matrix(total_effects[[1]])

    if (is_tvp) {
      # TVP model: use print_dynamics_tvp
      cat("\n")
      cat(strrep("=", 70), "\n")
      cat("FIT OBJECT: Time-Varying Parameter Model\n")
      cat(strrep("=", 70), "\n")

      # Print for each subject
      for (i in seq_along(total_effects)) {
        if (length(total_effects) > 1) {
          cat("\n")
          cat(strrep("-", 70), "\n")
          cat("Subject ", i, "\n")
          cat(strrep("-", 70), "\n")
        }
        print_dynamics_tvp(total_effects[[i]], period_labels = period_labels,
                          digits = digits, zero_char = zero_char,
                          highlight_nonzero = highlight_nonzero,
                          max_print = max_print, time_labels = time_labels)
      }

      return(invisible(fit_obj))
    } else {
      # Non-TVP model: extract total effects
      cat("\n")
      cat(strrep("=", 70), "\n")
      cat("FIT OBJECT: Total Effects\n")
      cat(strrep("=", 70), "\n")

      dynamics <- total_effects
      # Continue to regular processing below
    }
  }

  # Helper to print a single matrix
  print_single_matrix <- function(mat, label = NULL) {
    if (!is.matrix(mat) || !is.numeric(mat)) {
      stop("Input must be a numeric matrix")
    }

    d <- nrow(mat)
    p <- ncol(mat)

    # Print header
    if (!is.null(label)) {
      cat("\n")
      cat(strrep("=", 70), "\n")
      cat(label, "\n")
      cat(strrep("=", 70), "\n\n")
    }

    # Create row and column names if not present
    if (is.null(rownames(mat))) {
      if (time_labels) {
        rownames(mat) <- paste0("V", 1:d, "_t")
      } else {
        rownames(mat) <- paste0("V", 1:d)
      }
    }

    if (is.null(colnames(mat))) {
      if (time_labels) {
        colnames(mat) <- paste0("V", 1:p, "_{t-1}")
      } else {
        colnames(mat) <- paste0("V", 1:p)
      }
    }

    # Check if matrix is too large
    truncate_rows <- d > max_print
    truncate_cols <- p > max_print

    # Determine what to print
    if (truncate_rows || truncate_cols) {
      # For large matrices, print summary
      print_large_matrix_summary(mat, digits, zero_char, highlight_nonzero)
    } else {
      # For small matrices, print full
      print_full_matrix(mat, digits, zero_char, highlight_nonzero)
    }

    cat("\n")
  }

  # Helper to print full matrix
  print_full_matrix <- function(mat, digits, zero_char, highlight_nonzero) {
    d <- nrow(mat)
    p <- ncol(mat)

    # Format values
    mat_formatted <- matrix("", nrow = d, ncol = p)
    for (i in 1:d) {
      for (j in 1:p) {
        val <- mat[i, j]
        if (abs(val) < 1e-10) {
          mat_formatted[i, j] <- zero_char
        } else {
          formatted_val <- sprintf(paste0("%.", digits, "f"), val)
          # Add space before positive numbers for alignment with negatives
          if (val >= 0) {
            formatted_val <- paste0(" ", formatted_val)
          }
          mat_formatted[i, j] <- formatted_val
        }
      }
    }

    # Add row and column names
    rownames(mat_formatted) <- rownames(mat)
    colnames(mat_formatted) <- colnames(mat)

    # Print using standard print
    print(noquote(mat_formatted))

    # Print summary statistics
    cat("\n")
    num_nonzero <- sum(abs(mat) > 1e-10)
    num_total <- d * p
    sparsity <- 1 - (num_nonzero / num_total)
    cat(sprintf("Non-zero entries: %d / %d (%.1f%% sparse)\n",
                num_nonzero, num_total, sparsity * 100))

    # Show range of non-zero values
    if (num_nonzero > 0) {
      nonzero_vals <- mat[abs(mat) > 1e-10]
      cat(sprintf("Range of non-zero values: [%.3f, %.3f]\n",
                  min(nonzero_vals), max(nonzero_vals)))
    }
  }

  # Helper to print summary for large matrices
  print_large_matrix_summary <- function(mat, digits, zero_char, highlight_nonzero) {
    d <- nrow(mat)
    p <- ncol(mat)

    cat(sprintf("Large matrix (%d x %d)\n\n", d, p))

    # Show top-left corner
    show_rows <- min(10, d)
    show_cols <- min(10, p)

    cat("Top-left corner:\n")
    corner <- mat[1:show_rows, 1:show_cols, drop = FALSE]
    print_full_matrix(corner, digits, zero_char, highlight_nonzero)

    if (d > show_rows || p > show_cols) {
      cat("\n... (matrix truncated, showing first ", show_rows, " rows and ", show_cols, " columns)\n", sep = "")
    }

    # Overall statistics
    cat("\n")
    num_nonzero <- sum(abs(mat) > 1e-10)
    num_total <- d * p
    sparsity <- 1 - (num_nonzero / num_total)
    cat(sprintf("Full matrix: %d x %d\n", d, p))
    cat(sprintf("Non-zero entries: %d / %d (%.1f%% sparse)\n",
                num_nonzero, num_total, sparsity * 100))

    if (num_nonzero > 0) {
      nonzero_vals <- mat[abs(mat) > 1e-10]
      cat(sprintf("Range of non-zero values: [%.3f, %.3f]\n",
                  min(nonzero_vals), max(nonzero_vals)))

      # Show top 5 largest absolute values
      top_vals <- sort(abs(nonzero_vals), decreasing = TRUE)[1:min(5, length(nonzero_vals))]
      cat("\nTop 5 largest absolute values:\n")
      for (i in seq_along(top_vals)) {
        # Find position of this value
        idx <- which(abs(mat) == top_vals[i], arr.ind = TRUE)[1, , drop = FALSE]
        cat(sprintf("  %d. %.3f at [%s, %s]\n",
                    i, mat[idx[1], idx[2]],
                    rownames(mat)[idx[1]], colnames(mat)[idx[2]]))
      }
    }
  }

  # Main logic: handle single matrix or list
  if (is.matrix(dynamics)) {
    # Single matrix
    label <- if (!is.null(subject_label)) subject_label else "Dynamics Matrix"
    print_single_matrix(dynamics, label = label)
  } else if (is.list(dynamics)) {
    # List of matrices
    for (i in seq_along(dynamics)) {
      mat <- dynamics[[i]]
      if (is.matrix(mat)) {
        label <- if (!is.null(subject_label)) {
          paste0(subject_label, " - Subject ", i)
        } else {
          paste0("Subject ", i, " Dynamics")
        }
        print_single_matrix(mat, label = label)
      } else {
        warning("Element ", i, " is not a matrix, skipping")
      }
    }
  } else {
    stop("Input must be a matrix or a list of matrices")
  }

  invisible(dynamics)
}


#' Print time-varying dynamics (TVP models)
#'
#' Prints dynamics for TVP models where dynamics vary across time periods.
#' This is a specialized wrapper around print_dynamics() for TVP structures.
#'
#' @param dynamics_tvp A list of matrices (one per period), or a list of lists
#'   (for multiple subjects, each containing a list of period-specific matrices)
#' @param period_labels Optional character vector of period labels (e.g., c("Pre", "During", "Post"))
#' @param ... Additional arguments passed to print_dynamics()
#'
#' @return Invisibly returns the input (for chaining)
#' @export
#'
#' @examples
#' # Three periods with different dynamics
#' phi_period1 <- matrix(c(0.3, 0, 0.4, 0, 0.4, 0.5, 0, 0.3, 0, 0.3, 0.2, 0, 0, 0, 0, 0.4), 4, 4)
#' phi_period2 <- matrix(c(0.4, 0.3, 0, 0.2, 0, 0.6, 0.4, 0, 0.3, 0, 0.3, 0, 0, 0.4, 0, 0.5), 4, 4)
#' phi_period3 <- matrix(c(0.5, 0.2, 0, 0, 0, 0.4, 0.5, 0, 0.5, 0, 0.4, 0.2, 0, 0.5, 0, 0.3), 4, 4)
#'
#' # For single subject
#' print_dynamics_tvp(list(phi_period1, phi_period2, phi_period3))
#'
#' # With custom period labels
#' print_dynamics_tvp(list(phi_period1, phi_period2, phi_period3),
#'                    period_labels = c("Baseline", "Treatment", "Follow-up"))
print_dynamics_tvp <- function(dynamics_tvp, period_labels = NULL, ...) {

  # Detect structure: single subject or multiple subjects
  if (is.list(dynamics_tvp) && length(dynamics_tvp) > 0) {
    # Check if first element is a matrix (single subject) or list (multiple subjects)
    if (is.matrix(dynamics_tvp[[1]])) {
      # Single subject: list of matrices
      print_tvp_single_subject(dynamics_tvp, period_labels, ...)
    } else if (is.list(dynamics_tvp[[1]])) {
      # Multiple subjects: list of lists
      print_tvp_multiple_subjects(dynamics_tvp, period_labels, ...)
    } else {
      stop("Unrecognized TVP structure")
    }
  } else {
    stop("dynamics_tvp must be a list")
  }

  invisible(dynamics_tvp)
}


# Helper for single subject TVP
print_tvp_single_subject <- function(period_list, period_labels, ...) {
  n_periods <- length(period_list)

  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("TIME-VARYING DYNAMICS (", n_periods, " periods)\n", sep = "")
  cat(strrep("=", 70), "\n")

  for (p in seq_along(period_list)) {
    mat <- period_list[[p]]

    if (is.matrix(mat)) {
      # Create label
      if (!is.null(period_labels) && length(period_labels) >= p) {
        label <- paste0("Period ", p, ": ", period_labels[p])
      } else {
        label <- paste0("Period ", p)
      }

      print_dynamics(mat, subject_label = label, ...)
    }
  }
}


# Helper for multiple subjects TVP
print_tvp_multiple_subjects <- function(subject_list, period_labels, ...) {
  n_subjects <- length(subject_list)

  for (s in seq_along(subject_list)) {
    period_list <- subject_list[[s]]

    if (is.list(period_list)) {
      n_periods <- length(period_list)

      cat("\n")
      cat(strrep("=", 70), "\n")
      cat("SUBJECT ", s, " TIME-VARYING DYNAMICS (", n_periods, " periods)\n", sep = "")
      cat(strrep("=", 70), "\n")

      for (p in seq_along(period_list)) {
        mat <- period_list[[p]]

        if (is.matrix(mat)) {
          # Create label
          if (!is.null(period_labels) && length(period_labels) >= p) {
            label <- paste0("Subject ", s, " - Period ", p, ": ", period_labels[p])
          } else {
            label <- paste0("Subject ", s, " - Period ", p)
          }

          print_dynamics(mat, subject_label = label, ...)
        }
      }
    }
  }
}


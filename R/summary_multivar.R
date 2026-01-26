#' Summary for multivar fit objects
#'
#' Prints a summary of a fitted multivar model, including dataset information,
#' model specifications, and key results.
#'
#' @param fit A fitted multivar object (output from cv.multivar())
#'
#' @return Invisibly returns the input object
#' @export
#'
#' @examples
#' # fit <- cv.multivar(object)
#' # summary_multivar(fit)
summary_multivar <- function(fit) {

  # Extract information from the stored model object
  obj <- fit$obj

  if (is.null(obj)) {
    stop("Fit object does not contain the original model specification (missing $obj)")
  }

  # Basic information
  k <- obj@k
  n <- obj@n
  d <- obj@d

  # Model specifications
  tvp <- obj@tvp
  intercept <- obj@intercept

  # Header
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("MULTIVAR Model Summary\n")
  cat(strrep("=", 70), "\n\n")

  # Dataset Information
  cat("Dataset Information:\n")
  cat(strrep("-", 70), "\n")

  cat(sprintf("  Number of datasets (subjects): %d\n", k))

  # Number of timepoints
  if (length(n) == 1) {
    cat(sprintf("  Number of timepoints:          %d\n", n))
  } else {
    # Check if all equal
    if (all(n == n[1])) {
      cat(sprintf("  Number of timepoints:          %d (all subjects)\n", n[1]))
    } else {
      cat(sprintf("  Number of timepoints:          %s (varies by subject)\n",
                  paste(n, collapse = ", ")))
      cat(sprintf("                                 Range: [%d, %d]\n",
                  min(n), max(n)))
      cat(sprintf("                                 Mean: %.1f\n", mean(n)))
    }
  }

  # Number of variables
  if (length(d) == 1) {
    cat(sprintf("  Number of variables:           %d\n", d))
  } else {
    # Check if all equal
    if (all(d == d[1])) {
      cat(sprintf("  Number of variables:           %d (all subjects)\n", d[1]))
    } else {
      cat(sprintf("  Number of variables:           %s (varies by subject)\n",
                  paste(d, collapse = ", ")))
    }
  }

  cat("\n")

  # Model Specification
  cat("Model Specification:\n")
  cat(strrep("-", 70), "\n")

  # TVP information
  if (tvp) {
    cat("  Time-varying parameters:       YES\n")

    # Break information
    if (!is.null(obj@breaks) && length(obj@breaks) > 0) {
      # Count number of periods
      n_periods <- length(obj@breaks[[1]])
      cat(sprintf("  Number of periods:             %d\n", n_periods))

      # Show period ranges for each subject
      if (k == 1) {
        # Single subject: show period ranges
        period_info <- sapply(obj@breaks[[1]], function(x) {
          sprintf("t=%d-%d", min(x), max(x))
        })
        cat(sprintf("  Period ranges:                 %s\n",
                    paste(period_info, collapse = ", ")))
      } else {
        # Multiple subjects: show ranges by subject
        cat("  Period ranges by subject:\n")
        for (i in seq_len(k)) {
          if (i <= length(obj@breaks)) {
            period_info <- sapply(obj@breaks[[i]], function(x) {
              sprintf("t=%d-%d", min(x), max(x))
            })
            cat(sprintf("    Subject %d:                    %s\n",
                        i, paste(period_info, collapse = ", ")))
          }
        }
      }
    }
  } else {
    cat("  Time-varying parameters:       NO (stationary)\n")
  }

  # Intercept information
  if (intercept) {
    cat("  Intercepts:                    YES\n")
    if (!is.null(obj@pen_common_intercept)) {
      cat(sprintf("    - Common intercept penalized: %s\n",
                  ifelse(obj@pen_common_intercept, "YES", "NO")))
    }
    if (!is.null(obj@pen_unique_intercept)) {
      cat(sprintf("    - Unique intercept penalized: %s\n",
                  ifelse(obj@pen_unique_intercept, "YES", "NO")))
    }
  } else {
    cat("  Intercepts:                    NO\n")
  }

  cat("\n")

  # Adaptive Weights
  cat("Adaptive Weights:\n")
  cat(strrep("-", 70), "\n")

  # Estimation method
  weightest <- obj@weightest
  weightest_label <- toupper(weightest)
  cat(sprintf("  Method:                        %s\n", weightest_label))

  # Adaptive or standard
  if (!is.null(obj@lassotype)) {
    lassotype_label <- ifelse(obj@lassotype == "adaptive", "Adaptive", "Standard")
    cat(sprintf("  Type:                          %s\n", lassotype_label))
  }

  # Lambda selection
  if (!is.null(obj@lambda_choice)) {
    lambda_choice <- obj@lambda_choice
    lambda_label <- ifelse(lambda_choice == "lambda.min",
                           "Minimum CV error",
                           "1 SE rule")
    cat(sprintf("  Lambda selection:              %s (%s)\n",
                lambda_label, lambda_choice))
  }

  # Cross-validation
  if (!is.null(obj@cv)) {
    cv_type <- obj@cv
    cv_label <- ifelse(cv_type == "blocked", "Blocked", "Standard")
    cat(sprintf("  Cross-validation:              %s", cv_label))
    if (!is.null(obj@nfolds)) {
      cat(sprintf(" (%d-fold)\n", obj@nfolds))
    } else {
      cat("\n")
    }
  }

  cat("\n")

  # Overview of Dynamics
  cat("Overview of Dynamics:\n")
  cat(strrep("-", 70), "\n")

  # Helper to compute matrix statistics
  compute_matrix_stats <- function(mat, eps = 1e-10) {
    if (!is.matrix(mat)) return(NULL)
    nz <- sum(abs(mat) > eps)
    total <- length(mat)
    sparsity <- (1 - nz / total) * 100
    list(
      nonzero = nz,
      total = total,
      sparsity = sparsity,
      range = if (nz > 0) range(mat[abs(mat) > eps]) else c(0, 0)
    )
  }

  # Get dynamics matrices
  total_mats <- fit$mats$total
  common_mat <- fit$mats$common

  if (tvp) {
    # TVP model: show per-period statistics
    if (k == 1) {
      # Single subject TVP
      period_list <- total_mats[[1]]
      n_periods <- length(period_list)

      cat(sprintf("  Total effects (time-varying, %d periods):\n", n_periods))

      for (p in seq_len(n_periods)) {
        stats <- compute_matrix_stats(period_list[[p]])
        if (!is.null(stats)) {
          cat(sprintf("    Period %d:  %2d edges, %.1f%% sparse, range [%.3f, %.3f]\n",
                      p, stats$nonzero, stats$sparsity,
                      stats$range[1], stats$range[2]))
        }
      }
    } else {
      # Multiple subjects TVP
      cat(sprintf("  Total effects (time-varying, K=%d subjects):\n", k))

      for (i in seq_len(k)) {
        period_list <- total_mats[[i]]
        n_periods <- length(period_list)

        cat(sprintf("    Subject %d:\n", i))
        for (p in seq_len(n_periods)) {
          stats <- compute_matrix_stats(period_list[[p]])
          if (!is.null(stats)) {
            cat(sprintf("      Period %d:  %2d edges, %.1f%% sparse\n",
                        p, stats$nonzero, stats$sparsity))
          }
        }
      }
    }
  } else {
    # Non-TVP model: show overall statistics
    if (k == 1) {
      # Single subject
      stats <- compute_matrix_stats(total_mats[[1]])
      if (!is.null(stats)) {
        cat(sprintf("  Total effects:  %2d edges, %.1f%% sparse, range [%.3f, %.3f]\n",
                    stats$nonzero, stats$sparsity,
                    stats$range[1], stats$range[2]))
      }
    } else {
      # Multiple subjects: show common and individual
      if (!is.null(common_mat)) {
        stats_common <- compute_matrix_stats(common_mat)
        if (!is.null(stats_common)) {
          cat(sprintf("  Common effects: %2d edges, %.1f%% sparse, range [%.3f, %.3f]\n",
                      stats_common$nonzero, stats_common$sparsity,
                      stats_common$range[1], stats_common$range[2]))
        }
      }

      cat("\n  Total effects by subject:\n")
      for (i in seq_len(k)) {
        stats <- compute_matrix_stats(total_mats[[i]])
        if (!is.null(stats)) {
          cat(sprintf("    Subject %d:  %2d edges, %.1f%% sparse, range [%.3f, %.3f]\n",
                      i, stats$nonzero, stats$sparsity,
                      stats$range[1], stats$range[2]))
        }
      }
    }
  }

  cat("\n")
  cat(strrep("=", 70), "\n\n")

  invisible(fit)
}

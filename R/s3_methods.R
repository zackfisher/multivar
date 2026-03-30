#' S3 Methods for multivar_fit Objects
#'
#' Standard S3 methods for fitted multivar models.
#'
#' @name multivar_fit_methods
NULL

#' @describeIn multivar_fit_methods Print a fitted multivar model.
#' @param x A multivar_fit object (output from cv.multivar()).
#' @param ... Additional arguments (currently unused).
#' @method print multivar_fit
#' @export
print.multivar_fit <- function(x, ...) {
  cat("\n")
  cat("Fitted multivar Model\n")
  cat(strrep("=", 70), "\n\n")

  obj <- x$obj

  # Basic info
  cat("Model Structure:\n")
  cat(sprintf("  Subjects (K):           %d\n", obj@k))
  cat(sprintf("  Variables (d):          %s\n",
              if(length(unique(obj@d)) == 1) obj@d[1] else paste(obj@d, collapse=", ")))
  cat(sprintf("  Timepoints (n):         %s\n",
              if(length(unique(obj@n)) == 1) obj@n[1] else paste(obj@n, collapse=", ")))
  cat(sprintf("  Time-varying (TVP):     %s\n", if(obj@tvp) "YES" else "NO"))
  cat(sprintf("  Intercepts:             %s\n", if(obj@intercept) "YES" else "NO"))

  # Penalty info
  cat("\nPenalization:\n")
  cat(sprintf("  Method:                 %s\n", toupper(obj@weightest)))
  cat(sprintf("  Type:                   %s\n",
              if(obj@lassotype == "adaptive") "Adaptive" else "Standard"))

  # Model fit
  cat("\nModel Fit:\n")
  sel <- if (!is.null(x$selection)) x$selection else "cv"
  if (sel == "ebic") {
    cat(sprintf("  Selection criterion:    EBIC (gamma=%.1f)\n", x$ebic$gamma))
    cat(sprintf("  EBIC value:             %.2f\n", x$ebic$values[x$ebic$best_idx]))
    best_B <- x$beta[, , x$ebic$best_idx]
    cat(sprintf("  Nonzero coefficients:   %d\n", sum(abs(best_B) > 0)))
  } else {
    cat(sprintf("  Selection criterion:    CV (MSFE)\n"))
    cat(sprintf("  Lambda selected:        %.6f\n", x$hyperparams$lambda))
    cat(sprintf("  Mean CV error:          %.6f\n", min(colMeans(x$MSFE))))
  }

  # Results
  if (obj@tvp) {
    n_periods <- length(x$mats$total[[1]])
    cat(sprintf("\nResults: %d periods, %d subjects\n", n_periods, obj@k))
  } else {
    cat(sprintf("\nResults: %d subjects\n", obj@k))
  }

  cat("\n")
  cat("Use summary() for detailed results\n")
  cat("Use plot() to visualize dynamics\n")
  cat("Use coef() to extract coefficients\n")
  cat("\n")

  invisible(x)
}

#' @describeIn multivar_fit_methods Summarize a fitted multivar model.
#' @param object A multivar_fit object (output from cv.multivar()).
#' @param ... Additional arguments (currently unused).
#' @method summary multivar_fit
#' @export
summary.multivar_fit <- function(object, ...) {
  # Call the existing summary_multivar function
  summary_multivar(object)
  invisible(object)
}

#' @describeIn multivar_fit_methods Plot fitted dynamics or prevalence.
#' @param x A multivar_fit object (output from cv.multivar()).
#' @param type Character; `"dynamics"` or `"prevalence"`.
#' @param ... Passed to the plotting helpers.
#' @method plot multivar_fit
#' @export
plot.multivar_fit <- function(x, type = c("dynamics", "prevalence"), ...) {
  type <- match.arg(type)

  if (type == "dynamics") {
    print_dynamics(x, ...)
  } else if (type == "prevalence") {
    if (x$obj@k > 1) {
      print_edge_prevalence(x, ...)
    } else {
      message("Edge prevalence plot requires K>1 subjects.")
      message("Showing dynamics instead.")
      print_dynamics(x, ...)
    }
  }

  invisible(x)
}

#' @describeIn multivar_fit_methods Extract coefficient matrices.
#' @param object A multivar_fit object (output from cv.multivar()).
#' @param type Character; one of `"total"`, `"common"`, `"unique"`.
#' @param subject Optional integer subject index.
#' @param period Optional integer period index (TVP only).
#' @param ... Additional arguments (unused).
#' @method coef multivar_fit
#' @export
coef.multivar_fit <- function(object, type = c("total", "common", "unique"),
                               subject = NULL, period = NULL, ...) {
  type <- match.arg(type)

  mats <- object$mats[[type]]

  obj <- object$obj

  # Handle TVP models
  if (obj@tvp) {
    if (is.null(period)) {
      # Return all periods
      if (is.null(subject)) {
        # All subjects, all periods
        return(mats)
      } else {
        # Specific subject, all periods
        if (subject > obj@k || subject < 1) {
          stop(sprintf("subject must be between 1 and %d", obj@k))
        }
        return(mats[[subject]])
      }
    } else {
      # Specific period
      n_periods <- length(mats[[1]])
      if (period > n_periods || period < 1) {
        stop(sprintf("period must be between 1 and %d", n_periods))
      }

      if (is.null(subject)) {
        # All subjects, specific period
        return(lapply(mats, function(x) x[[period]]))
      } else {
        # Specific subject, specific period
        if (subject > obj@k || subject < 1) {
          stop(sprintf("subject must be between 1 and %d", obj@k))
        }
        return(mats[[subject]][[period]])
      }
    }
  } else {
    # Non-TVP models
    if (type == "common") {
      # Common matrix is just one matrix
      return(mats)
    } else {
      # Total or unique - indexed by subject
      if (is.null(subject)) {
        # Return all subjects
        return(mats)
      } else {
        # Specific subject
        if (subject > obj@k || subject < 1) {
          stop(sprintf("subject must be between 1 and %d", obj@k))
        }
        return(mats[[subject]])
      }
    }
  }
}

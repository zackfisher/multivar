#' Helper Functions for Accessing multivar Results
#'
#' Simplified functions for extracting dynamics matrices from fitted multivar models.
#' These handle both TVP and non-TVP models transparently.
#'
#' @name get_dynamics_helpers
#' @param fit A fitted multivar object (output from cv.multivar())
#' @param subject Subject index (1 to K). If NULL, returns all subjects.
#' @param period Period index (for TVP models). If NULL, returns all periods.
#' @param type Type of effects: "total", "common", or "unique"
#'
#' @details
#' These helper functions simplify access to dynamics matrices by handling the
#' different structures of TVP and non-TVP models automatically.
#'
#' **For non-TVP models:**
#' - `fit$mats$common` is a single matrix (shared across all subjects)
#' - `fit$mats$unique` is a list of K matrices (one per subject)
#' - `fit$mats$total` is a list of K matrices (common + unique per subject)
#'
#' **For TVP models:**
#' - `fit$mats$common` is a single matrix (time-invariant, shared across all subjects)
#' - `fit$mats$tvp` is a list of K lists, each containing P period matrices (time-varying per subject)
#' - `fit$mats$total` is a list of K lists of P matrices (common + tvp per subject per period)
#'
#' The helper functions hide this complexity from the user.
#'
#' @examples
#' \dontrun{
#' # Non-TVP model
#' fit <- cv.multivar(object)
#' get_dynamics(fit, subject = 1)              # Subject 1 total effects
#' get_common_effects(fit)                     # Common effects
#' get_dynamics_all_subjects(fit)              # All subjects
#'
#' # TVP model
#' fit_tvp <- cv.multivar(object_tvp)
#' get_dynamics(fit_tvp, subject = 1, period = 2)    # Subject 1, Period 2
#' get_common_effects(fit_tvp, period = 2)           # Common for Period 2
#' get_dynamics_all_subjects(fit_tvp, period = 2)    # All subjects, Period 2
#' }
NULL

#' @rdname get_dynamics_helpers
#' @export
get_dynamics <- function(fit, subject = NULL, period = NULL, type = c("total", "common", "unique")) {
  type <- match.arg(type)

  # Check if fit object is valid
  if (!inherits(fit, "multivar_fit") && !is.list(fit)) {
    stop("fit must be a multivar_fit object from cv.multivar()")
  }

  if (is.null(fit$mats)) {
    stop("fit object does not contain $mats component")
  }

  obj <- fit$obj
  if (is.null(obj)) {
    stop("fit object does not contain $obj component")
  }

  is_tvp <- obj@tvp
  k <- obj@k

  # Get the appropriate matrix type
  mats <- fit$mats[[type]]

  if (is.null(mats)) {
    stop(sprintf("fit object does not contain %s effects", type))
  }

  # Handle different cases
  if (type == "common") {
    # Common effects
    if (is_tvp) {
      # TVP: Common effects are time-invariant (single matrix)
      # Time-varying components are in fit$mats$tvp (per-subject, per-period)
      # Total[subj][period] = common + tvp[subj][period]
      # For period-specific total effects, use type="total" with period parameter
      if (!is.null(period)) {
        message("Note: For TVP models, 'common' effects are time-invariant (same for all periods).\n",
                "      Use get_total_effects() with period parameter for period-specific dynamics.")
      }
      return(mats)  # Return time-invariant common matrix
    } else {
      # Non-TVP: common is a single matrix
      if (!is.null(period)) {
        warning("period parameter ignored for non-TVP models")
      }
      return(mats)
    }
  } else {
    # Total or unique effects (indexed by subject)
    if (is_tvp) {
      # TVP: mats is a list of K lists of period matrices
      if (is.null(subject)) {
        # Return all subjects
        if (is.null(period)) {
          return(mats)  # All subjects, all periods
        } else {
          # All subjects, specific period
          n_periods <- length(mats[[1]])
          if (period < 1 || period > n_periods) {
            stop(sprintf("period must be between 1 and %d", n_periods))
          }
          return(lapply(mats, function(x) x[[period]]))
        }
      } else {
        # Specific subject
        if (subject < 1 || subject > k) {
          stop(sprintf("subject must be between 1 and %d", k))
        }

        if (is.null(period)) {
          return(mats[[subject]])  # All periods for this subject
        } else {
          n_periods <- length(mats[[subject]])
          if (period < 1 || period > n_periods) {
            stop(sprintf("period must be between 1 and %d", n_periods))
          }
          return(mats[[subject]][[period]])
        }
      }
    } else {
      # Non-TVP: mats is a list of K matrices
      if (!is.null(period)) {
        warning("period parameter ignored for non-TVP models")
      }

      if (is.null(subject)) {
        return(mats)  # All subjects
      } else {
        if (subject < 1 || subject > k) {
          stop(sprintf("subject must be between 1 and %d", k))
        }
        return(mats[[subject]])
      }
    }
  }
}

#' @rdname get_dynamics_helpers
#' @export
get_common_effects <- function(fit, period = NULL) {
  get_dynamics(fit, subject = NULL, period = period, type = "common")
}

#' @rdname get_dynamics_helpers
#' @export
get_unique_effects <- function(fit, subject = NULL, period = NULL) {
  get_dynamics(fit, subject = subject, period = period, type = "unique")
}

#' @rdname get_dynamics_helpers
#' @export
get_total_effects <- function(fit, subject = NULL, period = NULL) {
  get_dynamics(fit, subject = subject, period = period, type = "total")
}

#' @rdname get_dynamics_helpers
#' @export
get_dynamics_all_subjects <- function(fit, period = NULL, type = c("total", "common", "unique")) {
  type <- match.arg(type)
  get_dynamics(fit, subject = NULL, period = period, type = type)
}

#' Get Number of Periods in TVP Model
#'
#' Returns the number of time-varying periods in a TVP model.
#' Returns 1 for non-TVP models.
#'
#' @param fit A fitted multivar object
#' @return Integer number of periods
#' @export
#'
#' @examples
#' \dontrun{
#' n_periods <- get_n_periods(fit)
#' }
get_n_periods <- function(fit) {
  if (!inherits(fit, "multivar_fit") && !is.list(fit)) {
    stop("fit must be a multivar_fit object from cv.multivar()")
  }

  obj <- fit$obj
  if (is.null(obj)) {
    stop("fit object does not contain $obj component")
  }

  if (!obj@tvp) {
    return(1)
  }

  # TVP model - get number of periods from total effects
  total <- fit$mats$total
  if (is.null(total) || length(total) == 0) {
    stop("Cannot determine number of periods")
  }

  # Total is a list of subjects, each containing a list of period matrices
  return(length(total[[1]]))
}

#' Check if Model is Time-Varying
#'
#' Returns TRUE if the model has time-varying parameters (TVP).
#'
#' @param fit A fitted multivar object
#' @return Logical indicating if model is TVP
#' @export
#'
#' @examples
#' \dontrun{
#' if (is_tvp(fit)) {
#'   cat("Model has time-varying parameters\n")
#' }
#' }
is_tvp <- function(fit) {
  if (!inherits(fit, "multivar_fit") && !is.list(fit)) {
    stop("fit must be a multivar_fit object from cv.multivar()")
  }

  obj <- fit$obj
  if (is.null(obj)) {
    stop("fit object does not contain $obj component")
  }

  return(obj@tvp)
}

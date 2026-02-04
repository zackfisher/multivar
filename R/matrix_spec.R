#' Build Design Matrix Specification
#'
#' Creates a specification object that defines the column and row structure
#' of the design matrix A. This is the single source of truth for all
#' downstream functions that need to know column indices.
#'
#' @param k Integer. Number of subjects.
#' @param d Integer. Number of variables.
#' @param n Integer vector. Timepoints per subject (length k).
#' @param tvp Logical. Whether time-varying parameters are used.
#' @param common_effects Logical. Whether common effects are included.
#' @param subgroup Logical. Whether subgroup structure is present.
#' @param common_tvp_effects Logical. Whether common TVP effects are included.
#' @param breaks List. Period definitions for TVP models (length k).
#' @param subgroup_membership Integer vector. Subgroup assignments (length k).
#'
#' @return A list with class "matrix_spec" containing:
#'   \itemize{
#'     \item{params: Model parameters (k, d, n, flags)}
#'     \item{cols: Column index ranges for each effect type}
#'     \item{rows: Row index ranges for subjects and periods}
#'   }
#'
#' @details
#' The design matrix A has the following column structure depending on model type:
#'
#' \strong{k=1, non-TVP:}
#' \preformatted{[d columns]}
#'
#' \strong{k=1, TVP, common_effects=TRUE:}
#' \preformatted{[d common | d*P unique_tvp]}
#'
#' \strong{k=1, TVP, common_effects=FALSE:}
#' \preformatted{[d*P unique_tvp]}
#'
#' \strong{k>1, non-TVP, no subgroup:}
#' \preformatted{[d common | d unique_1 | ... | d unique_k]}
#'
#' \strong{k>1, non-TVP, with subgroup:}
#' \preformatted{[d common | d subgrp_1 | ... | d subgrp_S | d unique_1 | ... | d unique_k]}
#'
#' \strong{k>1, TVP, common_effects=TRUE, common_tvp_effects=TRUE:}
#' \preformatted{[d common | d unique_1 | ... | d unique_k | d*P common_tvp | d*P unique_tvp_1 | ... | d*P unique_tvp_k]}
#'
#' @export
build_matrix_spec <- function(k, d, n, tvp = FALSE, common_effects = TRUE,
                               subgroup = FALSE, common_tvp_effects = FALSE,
                               breaks = NULL, subgroup_membership = NULL) {


  # Validate inputs
  stopifnot(k >= 1)
  stopifnot(d >= 1)
  stopifnot(length(n) == k)
  if (!common_effects && !tvp) {
    stop("common_effects = FALSE requires tvp = TRUE")
  }
  if (common_tvp_effects && !tvp) {
    stop("common_tvp_effects = TRUE requires tvp = TRUE")
  }
  if (common_tvp_effects && k == 1) {
    stop("common_tvp_effects = TRUE requires k > 1")
  }
  if (subgroup && is.null(subgroup_membership)) {
    stop("subgroup = TRUE requires subgroup_membership")
  }
  if (tvp && is.null(breaks)) {
    stop("tvp = TRUE requires breaks")
  }

  # Derive period structure from breaks
  # Note: n is already the number of rows per subject (ntk from constructModel)
  if (tvp) {
    num_periods <- sapply(breaks, length)
    period_lengths <- lapply(breaks, function(b) sapply(b, length))
  } else {
    num_periods <- rep(1L, k)
    period_lengths <- lapply(n, function(ni) ni)  # n is already row count
  }

  # Derive subgroup structure
  if (subgroup) {
    num_subgroups <- max(subgroup_membership)
  } else {
    num_subgroups <- 0L
    subgroup_membership <- NULL
  }

  # Build column specification
  cols <- build_column_spec(k, d, tvp, common_effects, subgroup,
                            common_tvp_effects, num_periods, num_subgroups)

  # Build row specification
  rows <- build_row_spec(k, n, tvp, breaks)

  # Assemble the spec

  spec <- list(
    params = list(
      k = k,
      d = d,
      n = n,
      tvp = tvp,
      common_effects = common_effects,
      subgroup = subgroup,
      common_tvp_effects = common_tvp_effects,
      num_periods = num_periods,
      period_lengths = period_lengths,
      num_subgroups = num_subgroups,
      subgroup_membership = subgroup_membership,
      breaks = breaks
    ),
    cols = cols,
    rows = rows
  )

  class(spec) <- c("matrix_spec", "list")
  spec
}


#' Build Column Specification
#'
#' Internal function to compute column index ranges.
#'
#' @keywords internal
build_column_spec <- function(k, d, tvp, common_effects, subgroup,
                               common_tvp_effects, num_periods, num_subgroups) {

  col_idx <- 1L  # Current column position

  # --- Common effects ---
  if (common_effects) {
    common <- c(col_idx, col_idx + d - 1L)
    col_idx <- col_idx + d
  } else {
    common <- NULL
  }

  # --- Subgroup effects ---
  if (subgroup && num_subgroups > 0) {
    subgrp <- vector("list", num_subgroups)
    for (s in seq_len(num_subgroups)) {
      subgrp[[s]] <- c(col_idx, col_idx + d - 1L)
      col_idx <- col_idx + d
    }
  } else {
    subgrp <- NULL
  }

  # --- Unique effects (per subject) ---
  # For k=1: no separate unique columns in the design matrix
  #   - k=1 non-TVP: common IS the total (unique points to common)
  #   - k=1 TVP: structure is [common | tvp_periods], no unique base block
  # For k>1: always have unique columns per subject
  if (k == 1) {
    if (!tvp && common_effects) {
      # k=1 non-TVP: common IS the total, unique points to same columns
      unique <- list(common)
    } else {
      # k=1 TVP: no unique base block (tvp_periods serve this role)
      unique <- NULL
    }
  } else {
    unique <- vector("list", k)
    for (i in seq_len(k)) {
      unique[[i]] <- c(col_idx, col_idx + d - 1L)
      col_idx <- col_idx + d
    }
  }

  # --- Common TVP effects ---
  if (tvp && common_tvp_effects) {
    # Assumes same number of periods for all subjects for common TVP
    P <- num_periods[1]
    common_tvp <- vector("list", P)
    for (p in seq_len(P)) {
      common_tvp[[p]] <- c(col_idx, col_idx + d - 1L)
      col_idx <- col_idx + d
    }
  } else {
    common_tvp <- NULL
  }

  # --- Unique TVP effects (per subject, per period) ---
  if (tvp) {
    unique_tvp <- vector("list", k)
    for (i in seq_len(k)) {
      P_i <- num_periods[i]
      unique_tvp[[i]] <- vector("list", P_i)
      for (p in seq_len(P_i)) {
        unique_tvp[[i]][[p]] <- c(col_idx, col_idx + d - 1L)
        col_idx <- col_idx + d
      }
    }
  } else {
    unique_tvp <- NULL
  }

  list(
    common = common,
    subgrp = subgrp,
    unique = unique,
    common_tvp = common_tvp,
    unique_tvp = unique_tvp,
    total = col_idx - 1L
  )
}


#' Build Row Specification
#'
#' Internal function to compute row index ranges.
#'
#' @param n Integer vector. Number of rows per subject (already accounts for lag).
#'          This should be ntk from constructModel, which is nrow(dat$b).
#'
#' @keywords internal
build_row_spec <- function(k, n, tvp, breaks) {

  # n is already the number of rows per subject (ntk = nrow(b) after lag applied)
  row_counts <- n

  # Subject row ranges
  subject <- vector("list", k)
  row_idx <- 1L
  for (i in seq_len(k)) {
    subject[[i]] <- c(row_idx, row_idx + row_counts[i] - 1L)
    row_idx <- row_idx + row_counts[i]
  }

  # Period row ranges (within each subject)
  if (tvp && !is.null(breaks)) {
    period <- vector("list", k)
    for (i in seq_len(k)) {
      P_i <- length(breaks[[i]])
      period[[i]] <- vector("list", P_i)

      # breaks[[i]] contains the actual row indices for each period
      # These are relative to the subject's data
      subject_start <- subject[[i]][1]
      for (p in seq_len(P_i)) {
        period_rows <- breaks[[i]][[p]]
        # Convert to absolute row indices
        period[[i]][[p]] <- c(
          subject_start + min(period_rows) - 1L,
          subject_start + max(period_rows) - 1L
        )
      }
    }
  } else {
    period <- NULL
  }

  list(
    subject = subject,
    period = period,
    total = row_idx - 1L
  )
}


# --- Accessor Functions ---

#' Get Common Effect Columns
#' @param spec A matrix_spec object
#' @return Integer vector c(start, end) or NULL
#' @export
get_common_cols <- function(spec) {

  spec$cols$common
}

#' Get Unique Effect Columns for a Subject
#' @param spec A matrix_spec object
#' @param subject Integer. Subject index (1 to k)
#' @return Integer vector c(start, end) or NULL (for k=1 TVP)
#' @export
get_unique_cols <- function(spec, subject) {
  if (is.null(spec$cols$unique)) return(NULL)
  spec$cols$unique[[subject]]
}

#' Get Subgroup Effect Columns
#' @param spec A matrix_spec object
#' @param subgroup Integer. Subgroup index (1 to S)
#' @return Integer vector c(start, end) or NULL
#' @export
get_subgrp_cols <- function(spec, subgroup) {
  if (is.null(spec$cols$subgrp)) return(NULL)
  spec$cols$subgrp[[subgroup]]
}

#' Get Subgroup Columns for a Subject
#' @param spec A matrix_spec object
#' @param subject Integer. Subject index (1 to k)
#' @return Integer vector c(start, end) or NULL
#' @export
get_subgrp_cols_for_subject <- function(spec, subject) {
  if (is.null(spec$cols$subgrp)) return(NULL)
  sg <- spec$params$subgroup_membership[subject]
  spec$cols$subgrp[[sg]]
}

#' Get Common TVP Columns for a Period
#' @param spec A matrix_spec object
#' @param period Integer. Period index (1 to P)
#' @return Integer vector c(start, end) or NULL
#' @export
get_common_tvp_cols <- function(spec, period) {
  if (is.null(spec$cols$common_tvp)) return(NULL)
  spec$cols$common_tvp[[period]]
}

#' Get Unique TVP Columns for a Subject and Period
#' @param spec A matrix_spec object
#' @param subject Integer. Subject index (1 to k)
#' @param period Integer. Period index (1 to P_i)
#' @return Integer vector c(start, end) or NULL
#' @export
get_unique_tvp_cols <- function(spec, subject, period) {
  if (is.null(spec$cols$unique_tvp)) return(NULL)
  spec$cols$unique_tvp[[subject]][[period]]
}

#' Get Row Range for a Subject
#' @param spec A matrix_spec object
#' @param subject Integer. Subject index (1 to k)
#' @return Integer vector c(start, end)
#' @export
get_subject_rows <- function(spec, subject) {
  spec$rows$subject[[subject]]
}

#' Get Row Range for a Subject's Period
#' @param spec A matrix_spec object
#' @param subject Integer. Subject index (1 to k)
#' @param period Integer. Period index (1 to P_i)
#' @return Integer vector c(start, end) or NULL
#' @export
get_period_rows <- function(spec, subject, period) {
  if (is.null(spec$rows$period)) return(NULL)
  spec$rows$period[[subject]][[period]]
}


# --- Utility Functions ---
#' Expand Column Range to Indices
#'
#' Converts a c(start, end) range to a full sequence.
#' @param range Integer vector c(start, end) or NULL
#' @return Integer vector or NULL
#' @export
expand_range <- function(range) {
  if (is.null(range)) return(NULL)
  seq(range[1], range[2])
}

#' Print method for matrix_spec
#' @param x A matrix_spec object
#' @param ... Additional arguments (ignored)
#' @export
print.matrix_spec <- function(x, ...) {
  cat("Matrix Specification\n")
  cat("====================\n\n")


  cat("Parameters:\n")
  cat(sprintf("  k = %d subjects\n", x$params$k))
  cat(sprintf("  d = %d variables\n", x$params$d))
  cat(sprintf("  n = %s timepoints\n", paste(x$params$n, collapse = ", ")))
  cat(sprintf("  tvp = %s\n", x$params$tvp))
  cat(sprintf("  common_effects = %s\n", x$params$common_effects))
  cat(sprintf("  subgroup = %s\n", x$params$subgroup))
  cat(sprintf("  common_tvp_effects = %s\n", x$params$common_tvp_effects))
  if (x$params$tvp) {
    cat(sprintf("  num_periods = %s\n", paste(x$params$num_periods, collapse = ", ")))
  }
  cat("\n")

  cat("Column Structure:\n")
  cat(sprintf("  Total columns: %d\n", x$cols$total))
  if (!is.null(x$cols$common)) {
    cat(sprintf("  Common: cols %d-%d\n", x$cols$common[1], x$cols$common[2]))
  }
  if (!is.null(x$cols$subgrp)) {
    for (s in seq_along(x$cols$subgrp)) {
      cat(sprintf("  Subgroup %d: cols %d-%d\n", s,
                  x$cols$subgrp[[s]][1], x$cols$subgrp[[s]][2]))
    }
  }
  if (!is.null(x$cols$unique)) {
    for (i in seq_along(x$cols$unique)) {
      cat(sprintf("  Unique[%d]: cols %d-%d\n", i,
                  x$cols$unique[[i]][1], x$cols$unique[[i]][2]))
    }
  }
  if (!is.null(x$cols$common_tvp)) {
    for (p in seq_along(x$cols$common_tvp)) {
      cat(sprintf("  Common TVP[%d]: cols %d-%d\n", p,
                  x$cols$common_tvp[[p]][1], x$cols$common_tvp[[p]][2]))
    }
  }
  if (!is.null(x$cols$unique_tvp)) {
    for (i in seq_along(x$cols$unique_tvp)) {
      for (p in seq_along(x$cols$unique_tvp[[i]])) {
        cat(sprintf("  Unique TVP[%d][%d]: cols %d-%d\n", i, p,
                    x$cols$unique_tvp[[i]][[p]][1], x$cols$unique_tvp[[i]][[p]][2]))
      }
    }
  }
  cat("\n")

  cat("Row Structure:\n")
  cat(sprintf("  Total rows: %d\n", x$rows$total))
  for (i in seq_along(x$rows$subject)) {
    cat(sprintf("  Subject %d: rows %d-%d\n", i,
                x$rows$subject[[i]][1], x$rows$subject[[i]][2]))
  }

  invisible(x)
}


#' Validate Matrix Specification
#'
#' Checks internal consistency of a matrix_spec object.
#'
#' @param spec A matrix_spec object
#' @return TRUE if valid, otherwise throws an error
#' @export
validate_matrix_spec <- function(spec) {


  # Check class

  if (!inherits(spec, "matrix_spec")) {
    stop("Object is not a matrix_spec")
  }

  # Check column indices don't overlap and are contiguous
  all_ranges <- list()
  if (!is.null(spec$cols$common)) {
    all_ranges <- c(all_ranges, list(spec$cols$common))
  }
  if (!is.null(spec$cols$subgrp)) {
    all_ranges <- c(all_ranges, spec$cols$subgrp)
  }
  # For k=1 non-TVP, unique points to common, so skip duplicate check
  # For k=1 TVP, unique is NULL
  if (!is.null(spec$cols$unique) &&
      !(spec$params$k == 1 && !spec$params$tvp && spec$params$common_effects)) {
    all_ranges <- c(all_ranges, spec$cols$unique)
  }
  if (!is.null(spec$cols$common_tvp)) {
    all_ranges <- c(all_ranges, spec$cols$common_tvp)
  }
  if (!is.null(spec$cols$unique_tvp)) {
    for (i in seq_along(spec$cols$unique_tvp)) {
      all_ranges <- c(all_ranges, spec$cols$unique_tvp[[i]])
    }
  }

  # Check contiguity
  if (length(all_ranges) > 0) {
    starts <- sapply(all_ranges, `[`, 1)
    ends <- sapply(all_ranges, `[`, 2)
    order_idx <- order(starts)
    starts <- starts[order_idx]
    ends <- ends[order_idx]

    if (starts[1] != 1) {
      stop("Column indices don't start at 1")
    }
    for (i in seq_along(starts)[-1]) {
      if (starts[i] != ends[i-1] + 1) {
        stop(sprintf("Gap in column indices between %d and %d", ends[i-1], starts[i]))
      }
    }
    if (ends[length(ends)] != spec$cols$total) {
      stop("Column indices don't match total")
    }
  }

  # Check row indices
  row_ends <- sapply(spec$rows$subject, `[`, 2)
  if (max(row_ends) != spec$rows$total) {
    stop("Row indices don't match total")
  }

  TRUE
}

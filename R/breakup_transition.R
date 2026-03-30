#' Extract Coefficient Matrices from Fitted Model
#'
#' Extracts common, unique, subgroup, and TVP coefficient matrices from the
#' fitted coefficient array B using the matrix specification for column indices.
#'
#' @param B Matrix. The fitted coefficient matrix (d x total_cols).
#' @param spec List. Matrix specification from build_matrix_spec().
#' @param Ak List. List of design matrices (for column/row names).
#' @param breaks List. Period definitions for TVP models.
#'
#' @return List containing:
#'   \itemize{
#'     \item{common: Common effects matrix (d x d)}
#'     \item{subgrp: List of subgroup effects matrices}
#'     \item{unique: List of unique effects matrices per subject}
#'     \item{tvp: List of time-varying effects per subject per timepoint}
#'     \item{common_tvp: List of common TVP effects per period}
#'     \item{total: List of total effects (sum of all components)}
#'   }
#'
#' @keywords internal
breakup_transition <- function(B, spec, Ak, breaks = NULL) {


  # Extract parameters from spec
  k <- spec$params$k
  d <- spec$params$d
  tvp <- spec$params$tvp
  common_effects <- spec$params$common_effects
  subgroup <- spec$params$subgroup
  common_tvp_effects <- spec$params$common_tvp_effects
  subgroup_membership <- spec$params$subgroup_membership
  ntk <- spec$params$n

  # Get variable names from Ak
  varnames <- colnames(Ak[[1]])
  if (is.null(varnames)) {
    varnames <- paste0("V", seq_len(d))
  }

  # Helper to set matrix names
  set_names <- function(mat) {
    colnames(mat) <- varnames
    rownames(mat) <- varnames
    mat
  }

  #---------------------------------------------------------------------------
  # k=1 case
  #---------------------------------------------------------------------------
  if (k == 1) {

    if (!tvp) {
      # k=1 non-TVP: B is the total effect matrix
      common_mat <- set_names(B)
      unique_mats <- list(B)
      total_mats <- list(set_names(B))
      subgrp_mats <- NULL
      tvp_mats <- NULL
      common_tvp_mats <- NULL

    } else {
      # k=1 TVP
      num_periods <- spec$params$num_periods[1]

      # Helper to expand period-level matrices to timepoint-level
      expand_to_timepoints <- function(period_mats) {
        lapply(seq_len(ntk[1]), function(t) {
          period_idx <- which(sapply(breaks[[1]], function(b) t %in% b))
          period_mats[[period_idx]]
        })
      }

      if (!common_effects) {
        # k=1 TVP, common_effects=FALSE: only TVP columns
        common_mat <- NULL
        unique_mats <- list(set_names(matrix(0, d, d)))

        # Extract TVP matrices using spec
        tvp_period_mats <- lapply(seq_len(num_periods), function(p) {
          cols <- expand_range(spec$cols$unique_tvp[[1]][[p]])
          set_names(B[, cols, drop = FALSE])
        })

        tvp_mats <- list(expand_to_timepoints(tvp_period_mats))
        total_mats <- tvp_mats
        common_tvp_mats <- NULL

      } else {
        # k=1 TVP, common_effects=TRUE: common + TVP columns
        # Extract common (base) matrix
        common_cols <- expand_range(spec$cols$common)
        common_mat <- set_names(B[, common_cols, drop = FALSE])

        # For k=1, unique is zero
        unique_mats <- list(set_names(matrix(0, d, d)))

        # Extract TVP matrices using spec
        tvp_period_mats <- lapply(seq_len(num_periods), function(p) {
          cols <- expand_range(spec$cols$unique_tvp[[1]][[p]])
          set_names(B[, cols, drop = FALSE])
        })

        tvp_mats <- list(expand_to_timepoints(tvp_period_mats))

        # Total = common + tvp for each timepoint
        total_mats <- list(
          lapply(seq_len(ntk[1]), function(t) {
            common_mat + tvp_mats[[1]][[t]]
          })
        )
        common_tvp_mats <- NULL
      }

      subgrp_mats <- NULL
    }

  #---------------------------------------------------------------------------
  # k>1 case
  #---------------------------------------------------------------------------
  } else {

    # Extract common matrix (if present)
    if (common_effects) {
      common_cols <- expand_range(spec$cols$common)
      common_mat <- set_names(B[, common_cols, drop = FALSE])
    } else {
      common_mat <- NULL
    }

    if (!subgroup && !tvp) {
      #-----------------------------------------------------------------------
      # k>1, non-TVP, no subgroup
      #-----------------------------------------------------------------------
      unique_mats <- lapply(seq_len(k), function(i) {
        cols <- expand_range(spec$cols$unique[[i]])
        set_names(B[, cols, drop = FALSE])
      })

      if (common_effects) {
        total_mats <- lapply(unique_mats, function(mat) mat + common_mat)
      } else {
        total_mats <- unique_mats
      }

      subgrp_mats <- NULL
      tvp_mats <- NULL
      common_tvp_mats <- NULL

    } else if (subgroup && !tvp) {
      #-----------------------------------------------------------------------
      # k>1, non-TVP, with subgroup
      #-----------------------------------------------------------------------
      # Extract subgroup matrices
      subgrp_mats <- lapply(seq_len(k), function(i) {
        cols <- expand_range(get_subgrp_cols_for_subject(spec, i))
        set_names(B[, cols, drop = FALSE])
      })

      # Extract unique matrices
      unique_mats <- lapply(seq_len(k), function(i) {
        cols <- expand_range(spec$cols$unique[[i]])
        set_names(B[, cols, drop = FALSE])
      })

      # Total = common + subgroup + unique
      total_mats <- lapply(seq_len(k), function(i) {
        unique_mats[[i]] + subgrp_mats[[i]] + common_mat
      })

      tvp_mats <- NULL
      common_tvp_mats <- NULL

    } else if (tvp) {
      #-----------------------------------------------------------------------
      # k>1, TVP
      #-----------------------------------------------------------------------
      num_periods <- spec$params$num_periods[1]

      # Extract unique base matrices
      unique_mats <- lapply(seq_len(k), function(i) {
        cols <- expand_range(spec$cols$unique[[i]])
        set_names(B[, cols, drop = FALSE])
      })

      # Extract common TVP matrices (if enabled)
      if (common_tvp_effects && !is.null(spec$cols$common_tvp)) {
        common_tvp_mats <- lapply(seq_len(num_periods), function(p) {
          cols <- expand_range(spec$cols$common_tvp[[p]])
          set_names(B[, cols, drop = FALSE])
        })
      } else {
        common_tvp_mats <- NULL
      }

      # Extract unique TVP matrices
      # Structure: list of k subjects, each containing list of timepoints
      unique_tvp_mats <- lapply(seq_len(k), function(i) {
        num_periods_i <- spec$params$num_periods[i]

        # Extract per-period matrices
        period_mats <- lapply(seq_len(num_periods_i), function(p) {
          cols <- expand_range(spec$cols$unique_tvp[[i]][[p]])
          B[, cols, drop = FALSE]
        })

        # Expand to timepoint level
        lapply(seq_len(ntk[i]), function(t) {
          period_idx <- which(sapply(breaks[[i]], function(b) t %in% b))
          set_names(period_mats[[period_idx]])
        })
      })

      tvp_mats <- unique_tvp_mats

      # Build total effects: base + TVP for each timepoint
      total_mats <- lapply(seq_len(k), function(i) {
        lapply(seq_len(ntk[i]), function(t) {
          # Determine period for timepoint t
          period_idx <- which(sapply(breaks[[i]], function(b) t %in% b))

          # Start with base effects
          if (common_effects) {
            total <- common_mat + unique_mats[[i]]
          } else {
            total <- unique_mats[[i]]
          }

          # Add common TVP (if enabled)
          if (common_tvp_effects && !is.null(common_tvp_mats)) {
            total <- total + common_tvp_mats[[period_idx]]
          }

          # Add unique TVP
          total + unique_tvp_mats[[i]][[t]]
        })
      })

      subgrp_mats <- NULL
    }
  }

  list(
    common = common_mat,
    subgrp = subgrp_mats,
    unique = unique_mats,
    tvp = tvp_mats,
    common_tvp = common_tvp_mats,
    total = total_mats
  )
}

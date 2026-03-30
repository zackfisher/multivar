#' Compute EBIC for an array of fitted coefficient matrices
#'
#' Core computation of the Extended Bayesian Information Criterion (EBIC)
#' for each scenario in a beta array. Used internally by \code{cv.multivar}
#' when \code{selection = "ebic"} and by \code{select_by_ebic} for post-hoc
#' reselection.
#'
#' @param beta_array 3D array (d x p x n_scenarios) of fitted coefficients.
#' @param Z Design matrix (transposed), p x n.
#' @param Y Response matrix (transposed), d x n.
#' @param d Number of response variables.
#' @param gamma EBIC tuning parameter. \code{gamma = 0} gives standard BIC;
#'   \code{gamma = 0.5} (default) is a moderate EBIC; \code{gamma = 1} is
#'   the most conservative (strongest sparsity preference).
#' @param spec Optional matrix_spec object. When provided, degrees of freedom
#'   are counted per coefficient position within each block (common, unique,
#'   TVP), following Chen & Chen (2008). Each d x d position that has any
#'   nonzero entry counts as one degree of freedom, regardless of how many
#'   subjects share that position. This avoids over-penalizing shared structure.
#'   When NULL, falls back to raw nonzero counting.
#'
#' @return Numeric vector of EBIC values, one per scenario.
#'
#' @details
#' With \code{spec = NULL} (legacy behavior), EBIC counts every nonzero entry
#' in B as one degree of freedom:
#' \deqn{\text{EBIC}_i = \sum_{j=1}^{d} n \log(\text{RSS}_{ij} / n) + k_i \log(n) + 2 \gamma \, k_i \log(p)}
#'
#' With a \code{spec} object, degrees of freedom are counted per position:
#' for each d x d coefficient position (i,j), the position contributes 1 df
#' if any block (common, any unique, any TVP) has a nonzero entry at that
#' position. The eBIC extension uses \code{p = d} (predictors per equation)
#' instead of the full design matrix width, and employs the Chen & Chen (2008)
#' \code{lchoose(p, df)} formulation per equation.
#'
#' @export
compute_ebic <- function(beta_array, Z, Y, d, gamma = 0.5, spec = NULL) {
  n <- ncol(Z)
  p_full <- nrow(Z)
  n_scenarios <- dim(beta_array)[3]
  ebic_vals <- numeric(n_scenarios)

  if (is.null(spec)) {
    # Legacy: raw nonzero counting
    for (i in seq_len(n_scenarios)) {
      B_i <- beta_array[, , i]
      resid <- Y - B_i %*% Z
      rss_per_eq <- pmax(rowSums(resid^2), .Machine$double.eps)
      k_i <- sum(abs(B_i) > 0)
      ebic_vals[i] <- sum(n * log(rss_per_eq / n)) + k_i * log(n) + 2 * gamma * k_i * log(p_full)
    }
  } else {
    # Grouped df: count nonzero positions per equation, using p = d
    for (i in seq_len(n_scenarios)) {
      B_i <- beta_array[, , i]
      resid <- Y - B_i %*% Z
      rss_per_eq <- pmax(rowSums(resid^2), .Machine$double.eps)

      # Count df per block using position-based grouping:
      # Each block (common, unique_1, ..., unique_k) contributes the number
      # of d x d positions with any nonzero entry. This gives credit for
      # shared structure (one common position serving K subjects = 1 df)
      # while still penalizing unique block activations independently.
      df_per_block <- count_block_df(B_i, d, spec)
      df_total <- sum(df_per_block)

      # eBIC extension: lchoose(p, df) per block, where p = d^2
      # (candidate positions per block). Sum across blocks.
      p_positions <- d * d
      ebic_extension <- sum(sapply(df_per_block, function(df_b) {
        lchoose(p_positions, min(df_b, p_positions))
      }))

      ebic_vals[i] <- sum(n * log(rss_per_eq / n)) +
        df_total * log(n) +
        2 * gamma * ebic_extension
    }
  }

  ebic_vals
}


#' Count degrees of freedom per block using position-based grouping
#'
#' For each column block (common, unique_1, ..., unique_k, subgroup, TVP),
#' counts the number of active d x d positions (those with any nonzero entry
#' across all d equations). Returns a vector with one df count per block.
#'
#' @param B Coefficient matrix (d x total_cols)
#' @param d Number of response variables (= number of predictors per position)
#' @param spec Matrix specification object from build_matrix_spec()
#'
#' @return Integer vector, one entry per block, giving the number of active
#'   positions (0..d^2) in that block.
#' @keywords internal
count_block_df <- function(B, d, spec) {

  k <- spec$params$k
  df_blocks <- integer(0)

  # Helper: count active positions in a set of columns
  count_active <- function(cols) {
    if (is.null(cols) || length(cols) == 0) return(0L)
    sum(colSums(abs(B[, cols, drop = FALSE])) > 0)
  }

  # Common block
  df_blocks <- c(df_blocks, count_active(get_common_cols(spec)))

  # Unique blocks
  for (s in seq_len(k)) {
    df_blocks <- c(df_blocks, count_active(get_unique_cols(spec, s)))
  }

  # Subgroup blocks
  if (spec$params$subgroup) {
    n_subgroups <- max(spec$params$subgroup_membership)
    for (s in seq_len(n_subgroups)) {
      df_blocks <- c(df_blocks, count_active(get_subgrp_cols(spec, s)))
    }
  }

  # TVP blocks
  if (spec$params$tvp) {
    # Common TVP
    if (!is.null(spec$cols$common_tvp)) {
      for (p_idx in seq_along(spec$cols$common_tvp)) {
        df_blocks <- c(df_blocks, count_active(spec$cols$common_tvp[[p_idx]]))
      }
    }
    # Unique TVP
    if (!is.null(spec$cols$unique_tvp)) {
      for (s in seq_len(k)) {
        for (p_idx in seq_along(spec$cols$unique_tvp[[s]])) {
          df_blocks <- c(df_blocks, count_active(spec$cols$unique_tvp[[s]][[p_idx]]))
        }
      }
    }
  }

  df_blocks
}


#' Reselect best model using (E)BIC instead of MSFE
#'
#' Takes the output of \code{cv.multivar()} and reselects the best penalty
#' scenario using the Extended Bayesian Information Criterion (EBIC) or
#' standard BIC. This can improve structure recovery (fewer false positives)
#' compared to prediction-based cross-validation, especially when n >> p.
#'
#' @param fit Output of \code{cv.multivar()}.
#' @param gamma EBIC tuning parameter. \code{gamma = 0} gives standard BIC;
#'   \code{gamma = 0.5} (default) is a moderate EBIC; \code{gamma = 1} is
#'   the most conservative (strongest sparsity preference).
#'
#' @return A modified copy of \code{fit} with:
#'   \item{mats}{Transition matrices from the EBIC-selected scenario}
#'   \item{ebic}{Named list: \code{values} (EBIC for all scenarios),
#'         \code{best_idx} (selected scenario index), \code{gamma} (tuning value)}
#'   The original MSFE results are preserved.
#'
#' @export
select_by_ebic <- function(fit, gamma = 0.5) {

  if (is.null(fit$beta)) {
    stop("select_by_ebic requires the full beta array. ",
         "Re-run cv.multivar() with save_beta = TRUE (the default).")
  }

  obj <- fit$obj

  # Reconstruct Z and Y (transposed, as used in cv_blocked)
  Z <- t(as.matrix(obj@A))
  Y <- t(as.matrix(obj@b))

  ebic_vals <- compute_ebic(fit$beta, Z, Y, obj@d, gamma, obj@spec)
  best_idx <- which.min(ebic_vals)

  # Rebuild mats from the EBIC-selected scenario
  mats <- breakup_transition(
    B = fit$beta[, , best_idx],
    spec = obj@spec,
    Ak = obj@Ak,
    breaks = obj@breaks
  )

  # Recover intercepts if needed
  if (obj@intercept && length(obj@data_means) > 0) {
    mats$intercepts <- recover_intercepts(
      mats = mats,
      data_means = obj@data_means,
      k = obj@k,
      d = obj@d,
      subgroup = obj@subgroup,
      subgroup_membership = obj@subgroup_membership,
      tvp = obj@tvp,
      breaks = obj@breaks
    )
  }

  # Return modified fit
  fit$mats <- mats
  fit$ebic <- list(
    values   = ebic_vals,
    best_idx = best_idx,
    gamma    = gamma
  )

  fit
}

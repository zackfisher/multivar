#' Recover Intercepts from Centered Data
#'
#' After fitting the model on centered data, this function recovers the intercepts
#' using the formula: c = mean(b) - Phi * mean(A). For multi-group models, it also
#' decomposes intercepts into common and group-specific components.
#'
#' @param mats List. The fitted coefficient matrices from breakup_transition, containing:
#'   - common: Common effects matrix (Psi)
#'   - unique: List of unique effects matrices (Psi^(k))
#'   - total: List of total effects matrices (Phi^(k) = Psi + Psi^(k))
#'   - subgrp: List of subgroup effects matrices (if subgroups)
#' @param data_means List. For each group, a list containing mean_b and mean_A from original data
#' @param k Numeric. Number of groups
#' @param d Numeric. Number of variables
#' @param subgroup Logical. Whether subgroup structure is present
#' @param subgroup_membership Numeric vector. Subgroup assignments for each subject
#' @param tvp Logical. Whether time-varying parameters are present
#' @param breaks List. Time period definitions for TVP models
#'
#' @return List containing:
#'   - intercepts_total: List of total intercepts for each group (c_total^(k))
#'   - intercept_common: Common intercept (c)
#'   - intercepts_unique: List of group-specific intercepts (c^(k))
#'   - intercepts_subgrp: List of subgroup-specific intercepts (if applicable)
#'   - intercepts_tvp: List of time-varying intercepts for TVP models (if applicable)
#'
#' @details
#' For each group k, the total intercept is recovered using:
#'   c_total^(k) = mean(b^(k)) - Phi^(k) * mean(A^(k))
#' where b^(k) are the outcomes (Y_t) and A^(k) are the predictors (Y_{t-1}).
#' This is the same approach used by glmnet for intercept recovery.
#'
#' For multi-group models (k > 1), intercepts are decomposed as:
#'   c = (1/K) * sum_k(c_total^(k))     # Common intercept
#'   c^(k) = c_total^(k) - c             # Group-specific intercept
#'
#' This decomposition ensures the constraint: sum_k(c^(k)) = 0
#'
#' Note: This formula c = mean(b) - Phi * mean(A) is different from the steady-state
#' formula c = (I - Phi) * mean(Y), which assumes mean(b) = mean(A) = mean(Y).
#' In finite samples they give similar but not identical results. The formula used
#' here is preferred as it matches standard practice in penalized regression (e.g., glmnet).
#'
#' For TVP models, this function recovers time-varying intercepts c_t for each period.
#' When data_means contains per-period means (mean_b_t, mean_A_t), time-specific intercepts
#' are computed. Otherwise, falls back to base intercepts.
#'
#' @export
recover_intercepts <- function(mats, data_means, k, d, subgroup = FALSE,
                                subgroup_membership = NULL, tvp = FALSE, breaks = NULL) {

  # If no data_means provided (intercept=FALSE), return NULL
  if (is.null(data_means) || length(data_means) == 0) {
    return(NULL)
  }

  # Identity matrix
  I_mat <- diag(d[1])

  # Handle TVP intercepts if applicable
  if (tvp && !is.null(breaks) && !is.null(mats$total)) {
    # For TVP models, compute time-varying intercepts c_t for each period
    # Check if data_means contains per-period means
    has_period_means <- !is.null(data_means[[1]]$mean_b_periods)

    if (has_period_means) {
      # Compute time-varying intercepts using per-period means
      intercepts_tvp <- lapply(seq_len(k), function(i) {
        n_periods <- length(breaks[[i]])
        period_intercepts <- vector("list", n_periods)

        for (t in seq_len(n_periods)) {
          # Get Phi_t for this period
          Phi_t <- mats$total[[i]][[t]]

          # Get means for this period
          mean_b_t <- data_means[[i]]$mean_b_periods[[t]]
          mean_A_t <- data_means[[i]]$mean_A_periods[[t]]

          # Transform back to original scale if standardized
          # Note: Both A and b are standardized using sd from A (see setup_data.R)
          # So transformation uses sd_A for both: diag(sd_A) * Phi * diag(1/sd_A) = Phi
          # No transformation needed when same scaling, but we keep structure for clarity
          if (!is.null(data_means[[i]]$sd_A_periods)) {
            sd_A_t <- data_means[[i]]$sd_A_periods[[t]]
            # Since both A and b scaled by same sd_A, this is identity: Phi_orig = Phi_std
            Phi_t <- diag(sd_A_t) %*% Phi_t %*% diag(1/sd_A_t)
          }

          # Compute c_t = mean(b_t) - Phi_t * mean(A_t)
          c_t <- mean_b_t - Phi_t %*% mean_A_t
          period_intercepts[[t]] <- as.vector(c_t)
        }

        period_intercepts
      })

      # For TVP, don't compute common/unique decomposition on intercepts
      # (could be added later if needed)
      return(list(
        intercepts_total = NULL,
        intercept_common = NULL,
        intercepts_unique = NULL,
        intercepts_subgrp = NULL,
        intercepts_tvp = intercepts_tvp
      ))
    }
    # If no period means, fall through to compute base intercepts
  }

  # Extract total effects matrices (non-TVP or TVP fallback)
  total_effects <- mats$total

  # If total effects not available, construct from common + unique
  if (is.null(total_effects) && !is.null(mats$common) && !is.null(mats$unique)) {
    if (subgroup && !is.null(mats$subgrp)) {
      # With subgroups: Phi^(k) = common + subgrp + unique
      total_effects <- lapply(seq_along(mats$unique), function(i) {
        mats$common + mats$subgrp[[i]] + mats$unique[[i]]
      })
    } else {
      # No subgroups: Phi^(k) = common + unique
      total_effects <- lapply(mats$unique, function(mat) {
        mats$common + mat
      })
    }
  }

  # For k=1, total effects might be stored directly
  if (k == 1 && is.matrix(total_effects)) {
    total_effects <- list(total_effects)
  }

  # Recover total intercept for each group using: c = mean(b) - Phi * mean(A)
  # This is the same formula glmnet uses for intercept recovery
  intercepts_total <- lapply(seq_len(k), function(i) {
    Phi_k <- total_effects[[i]]
    mean_b <- data_means[[i]]$mean_b
    mean_A <- data_means[[i]]$mean_A

    # If data was standardized, Phi_k is in standardized units
    # Transform it back to original scale: Phi_orig = diag(sd_b) * Phi_std * diag(1/sd_A)
    if (!is.null(data_means[[i]]$sd_b) && !is.null(data_means[[i]]$sd_A)) {
      sd_b <- data_means[[i]]$sd_b
      sd_A <- data_means[[i]]$sd_A
      Phi_k <- diag(sd_b) %*% Phi_k %*% diag(1/sd_A)
    }

    c_total <- mean_b - Phi_k %*% mean_A
    as.vector(c_total)
  })

  # For k=1, no decomposition needed
  if (k == 1) {
    return(list(
      intercepts_total = intercepts_total,
      intercept_common = NULL,
      intercepts_unique = NULL,
      intercepts_subgrp = NULL
    ))
  }

  # For k>1, decompose into common and group-specific
  # Common intercept: c = (1/K) * sum_k(c_total^(k))
  intercept_common <- Reduce("+", intercepts_total) / k

  # Group-specific intercepts: c^(k) = c_total^(k) - c
  intercepts_unique <- lapply(intercepts_total, function(c_total) {
    c_total - intercept_common
  })

  # Subgroup-specific intercepts (if applicable)
  intercepts_subgrp <- NULL
  if (subgroup && !is.null(mats$subgrp) && !is.null(subgroup_membership)) {
    # Decompose: c_total^(k) = c + c_subgrp^(s) + c_unique^(k)
    # where s = subgroup_membership[k]

    # Get unique subgroup IDs
    unique_subgroups <- unique(subgroup_membership)
    n_subgroups <- length(unique_subgroups)

    # Compute subgroup-specific intercepts by averaging within each subgroup
    subgrp_intercepts_list <- lapply(unique_subgroups, function(sg) {
      # Find subjects in this subgroup
      subjects_in_sg <- which(subgroup_membership == sg)

      # Average total intercepts for subjects in this subgroup
      sg_total_mean <- Reduce("+", intercepts_total[subjects_in_sg]) / length(subjects_in_sg)

      # Subgroup intercept = (mean of subgroup totals) - common
      sg_total_mean - intercept_common
    })

    # Replicate subgroup intercepts for each subject based on membership
    intercepts_subgrp <- lapply(seq_len(k), function(i) {
      sg_id <- subgroup_membership[i]
      sg_index <- which(unique_subgroups == sg_id)
      subgrp_intercepts_list[[sg_index]]
    })

    # Recompute unique intercepts accounting for subgroup effects
    # c_unique^(k) = c_total^(k) - c - c_subgrp^(s)
    intercepts_unique <- lapply(seq_len(k), function(i) {
      intercepts_total[[i]] - intercept_common - intercepts_subgrp[[i]]
    })
  }

  return(list(
    intercepts_total = intercepts_total,
    intercept_common = intercept_common,
    intercepts_unique = intercepts_unique,
    intercepts_subgrp = intercepts_subgrp
  ))
}

#' Compute adaptive LASSO weights from initial estimates
#'
#' @param coefs Matrix of initial coefficient estimates
#' @param adapower Power for adaptive weighting (default 1)
#' @param pendiag Whether to penalize diagonal (AR effects). If FALSE, diagonal gets weight 1e-10
#' @return Weight matrix with same dimensions as coefs
#' @keywords internal
adaptive_weights <- function(coefs, adapower = 1, pendiag = TRUE) {
  w <- 1/abs(coefs)^adapower
  w[is.infinite(w)] <- 1e10
  if (!pendiag) diag(w) <- 1e-10
  w
}

#' Replicate and scale weight matrix by ratio parameters
#'
#' Creates a 3D weight array by replicating the base weight matrix across
#' all ratio combinations, then scales appropriate column blocks by their
#' corresponding ratio values.
#'
#' @param w_mat Base weight matrix (d x n_cols)
#' @param k Number of subjects
#' @param n_pred Number of predictors per subject
#' @param subgroup Whether subgrouping is enabled
#' @param tvp Whether time-varying parameters are enabled
#' @param include_common_effects Whether common effects are included
#' @param common_tvp_effects Whether common TVP effects are included
#' @param ratios_unique Ratio values for unique effects
#' @param ratios_subgroup Ratio values for subgroup effects
#' @param ratios_unique_tvp Ratio values for unique TVP effects
#' @param ratios_common_tvp Ratio values for common TVP effects
#' @param subgroup_membership Vector of subgroup assignments
#' @param initcoefs Initial coefficients (for checking common_tvp_effects)
#'
#' @return 3D weight array (d x n_cols x n_scenarios)
#' @keywords internal
scale_weights_by_ratios <- function(w_mat, k, n_pred, subgroup, tvp,
                                     include_common_effects, common_tvp_effects,
                                     ratios_unique, ratios_subgroup,
                                     ratios_unique_tvp, ratios_common_tvp,
                                     subgroup_membership, initcoefs) {

 # Step 1: Determine number of scenarios and replicate w_mat
 if (tvp) {
   if (k == 1) {
     n_scenarios <- length(ratios_unique_tvp)
   } else if (include_common_effects && common_tvp_effects) {
     if (is.null(ratios_common_tvp)) ratios_common_tvp <- ratios_unique_tvp
     n_scenarios <- length(ratios_unique) * length(ratios_unique_tvp) * length(ratios_common_tvp)
   } else if (include_common_effects && !common_tvp_effects) {
     n_scenarios <- length(ratios_unique) * length(ratios_unique_tvp)
   } else if (!include_common_effects && common_tvp_effects) {
     if (is.null(ratios_common_tvp)) ratios_common_tvp <- ratios_unique_tvp
     n_scenarios <- length(ratios_unique_tvp) * length(ratios_common_tvp)
   } else {
     n_scenarios <- length(ratios_unique_tvp)
   }
 } else if (subgroup) {
   n_scenarios <- length(ratios_unique) * length(ratios_subgroup)
 } else {
   n_scenarios <- length(ratios_unique)
 }

 W <- replicate(n_scenarios, w_mat, simplify = "array")

 # Step 2: Scale columns by ratios
 # Note: assumes all subjects have same number of predictors

 if (k == 1) {
   # Single subject cases
   if (!tvp) {
     # Non-TVP k=1: scale all by ratios_unique
     for (r in seq_along(ratios_unique)) {
       W[, , r] <- W[, , r] * ratios_unique[r]
     }
   } else if (include_common_effects) {
     # TVP k=1 with common: scale common columns by ratios_unique_tvp
     for (a in seq_along(ratios_unique_tvp)) {
       W[, 1:n_pred, a] <- W[, 1:n_pred, a] * ratios_unique_tvp[a]
     }
   }
   # TVP k=1 without common: no scaling needed

 } else {
   # Multiple subject cases
   if (!subgroup && !tvp) {
     # Standard k>1: scale common columns by ratios_unique
     if (include_common_effects) {
       for (r in seq_along(ratios_unique)) {
         W[, 1:n_pred, r] <- W[, 1:n_pred, r] * ratios_unique[r]
       }
     }

   } else if (tvp) {
     cnt <- 1

     if (include_common_effects && common_tvp_effects) {
       # 4-layer: common base + unique base + common TVP + unique TVP
       if (is.null(ratios_common_tvp)) ratios_common_tvp <- ratios_unique_tvp

       common_base_cols <- 1:n_pred
       tvp_start <- n_pred * (k + 1) + 1
       total_tvp_cols <- ncol(w_mat) - n_pred * (k + 1)

       # Estimate column split for common vs unique TVP
       if (!is.null(initcoefs$common_tvp_effects)) {
         common_tvp_cols <- tvp_start:(tvp_start + total_tvp_cols %/% 3 - 1)
         unique_tvp_cols <- (max(common_tvp_cols) + 1):ncol(w_mat)
       } else {
         unique_tvp_cols <- tvp_start:ncol(w_mat)
         common_tvp_cols <- NULL
       }

       for (r in seq_along(ratios_unique)) {
         for (a in seq_along(ratios_unique_tvp)) {
           for (b in seq_along(ratios_common_tvp)) {
             W[, common_base_cols, cnt] <- W[, common_base_cols, cnt] * ratios_unique[r]
             if (!is.null(common_tvp_cols) && length(common_tvp_cols) > 0) {
               W[, common_tvp_cols, cnt] <- W[, common_tvp_cols, cnt] * ratios_common_tvp[b]
             }
             W[, unique_tvp_cols, cnt] <- W[, unique_tvp_cols, cnt] * ratios_unique_tvp[a]
             cnt <- cnt + 1
           }
         }
       }

     } else if (include_common_effects && !common_tvp_effects) {
       # 3-layer: common base + unique base + unique TVP
       tvp_cols <- (n_pred * (k + 1) + 1):ncol(W)

       for (r in seq_along(ratios_unique)) {
         for (j in seq_along(ratios_unique_tvp)) {
           W[, 1:n_pred, cnt] <- W[, 1:n_pred, cnt] * ratios_unique[r]
           W[, tvp_cols, cnt] <- W[, tvp_cols, cnt] * ratios_unique_tvp[j]
           cnt <- cnt + 1
         }
       }

     } else if (!include_common_effects && common_tvp_effects) {
       # 3-layer: unique base + common TVP + unique TVP
       if (is.null(ratios_common_tvp)) ratios_common_tvp <- ratios_unique_tvp

       tvp_start <- n_pred * k + 1
       total_tvp_cols <- ncol(w_mat) - n_pred * k
       common_tvp_cols <- tvp_start:(tvp_start + total_tvp_cols %/% 3 - 1)
       unique_tvp_cols <- (max(common_tvp_cols) + 1):ncol(w_mat)

       for (a in seq_along(ratios_unique_tvp)) {
         for (b in seq_along(ratios_common_tvp)) {
           W[, common_tvp_cols, cnt] <- W[, common_tvp_cols, cnt] * ratios_common_tvp[b]
           W[, unique_tvp_cols, cnt] <- W[, unique_tvp_cols, cnt] * ratios_unique_tvp[a]
           cnt <- cnt + 1
         }
       }
     }
     # 2-layer (!include_common_effects && !common_tvp_effects): no scaling

   } else if (subgroup && !tvp) {
     # Subgrouping: scale common by ratios_unique, subgroup by ratios_subgroup
     subgroup_end <- n_pred * max(subgroup_membership) + n_pred
     cnt <- 1

     for (r in seq_along(ratios_unique)) {
       for (j in seq_along(ratios_subgroup)) {
         W[, 1:n_pred, cnt] <- W[, 1:n_pred, cnt] * ratios_unique[r]
         W[, (n_pred + 1):subgroup_end, cnt] <- W[, (n_pred + 1):subgroup_end, cnt] * ratios_subgroup[j]
         cnt <- cnt + 1
       }
     }
   }
 }

 W
}

#' @export
est_base_weight_mat <- function(
  W,
  Ak,
  initcoefs,
  ratios_unique,
  d,
  k,
  lassotype,
  weightest,
  subgroup_membership,
  subgroup,
  ratios_subgroup,
  pendiag,
  tvp,
  ratios_unique_tvp,
  ratios_common_tvp = NULL,
  intercept,
  common_effects = TRUE,
  common_tvp_effects = TRUE,
  spec = NULL){

  # Note: if intercepts are included, pendiag should account for AR effects

  # Save parameter value before it could be overwritten

  include_common_effects <- common_effects


  n_pred_per_subj <- vapply(Ak, function(x) { ncol(x) }, numeric(1))

  # n_responses = number of response equations (always d)
  # n_predictors = number of predictors per equation
  # (These are equal now, but kept separate for future extensibility
  # if additional predictors are added to the model)
  n_responses <- d
  n_predictors <- n_pred_per_subj[1]

  adapower <- 1

  if (lassotype == "standard"){
   
    w_mat <- W
   
  } else {

    if(length(Ak) == 1){

      if(!tvp){
        # Non-TVP k=1 case
        w_mat <- adaptive_weights(initcoefs$total_effects[[1]], adapower, pendiag)

      } else {
        # TVP k=1 case

        if (!common_effects || is.null(initcoefs$common_effects)) {
          # No time-invariant component: weights from per-period estimates directly
          # tvp_effects[[1]] contains the per-period estimates (not deviations)
          tvp_list <- lapply(seq_along(initcoefs$tvp_effects[[1]]), function(j){
            adaptive_weights(initcoefs$tvp_effects[[1]][[j]], adapower, pendiag)
          })

          # Reformat tvp weights to match A matrix structure (block-diagonal by period)
          # Each row j contains: [period1_var1...period1_vard | period2_var1...period2_vard | ...]
          w_mat <- do.call(rbind, lapply(1:n_responses, function(j){
            unlist(lapply(seq_along(tvp_list), function(m) tvp_list[[m]][j,]))
          }))

        } else {
          # Time-invariant (common) weights from common_effects
          # (For k=1 TVP, common_effects = median of per-period estimates)
          common_weights <- adaptive_weights(initcoefs$common_effects, adapower, pendiag)

          # TVP weights from tvp_effects (period-specific deviations from common)
          tvp_list <- lapply(seq_along(initcoefs$tvp_effects[[1]]), function(j){
            adaptive_weights(initcoefs$tvp_effects[[1]][[j]], adapower, pendiag = TRUE)
          })

          # Reformat tvp weights to match A matrix structure (block-diagonal by period)
          # Each row j contains: [period1_var1...period1_vard | period2_var1...period2_vard | ...]
          phi_weights <- do.call(rbind, lapply(1:n_responses, function(j){
            unlist(lapply(seq_along(tvp_list), function(m) tvp_list[[m]][j,]))
          }))

          # Combine common and TVP weights
          w_mat <- cbind(common_weights, phi_weights)
        }
      }

    } else {
      
      if(!subgroup){

        unique_weights_list <- lapply(seq_along(Ak), function(i){
          adaptive_weights(initcoefs$unique_effects[[i]], adapower, pendiag = TRUE)
        })

        if (include_common_effects){
          # Standard: common + unique weights
          common_weights <- adaptive_weights(initcoefs$common_effects, adapower, pendiag)
          w_mat <- cbind(common_weights, do.call("cbind", unique_weights_list))
        } else {
          # No common effects: only unique weights
          w_mat <- do.call("cbind", unique_weights_list)
        }

      } else {

        subgroup_weights_list <- lapply(seq_along(initcoefs$subgroup_effects), function(i){
          adaptive_weights(initcoefs$subgroup_effects[[i]], adapower, pendiag = TRUE)
        })

        unique_weights_list <- lapply(seq_along(Ak), function(i){
          adaptive_weights(initcoefs$unique_effects[[i]], adapower, pendiag = TRUE)
        })

        common_weights <- adaptive_weights(initcoefs$common_effects, adapower, pendiag)
        w_mat <- cbind(common_weights, do.call("cbind", subgroup_weights_list), do.call("cbind", unique_weights_list))

      }
      
      # TODO: TVP-specific intercept handling not yet implemented
      if(tvp){

        # 1. Compute unique TVP weights
        tvp_list <- lapply(seq_along(initcoefs$tvp_effects), function(i){
          lapply(seq_along(initcoefs$tvp_effects[[i]]), function(j){
            adaptive_weights(initcoefs$tvp_effects[[i]][[j]], adapower, pendiag = TRUE)
          })
        })

        # 2. Compute unique base weights
        unique_weights_list <- lapply(seq_along(Ak), function(i){
          adaptive_weights(initcoefs$unique_effects[[i]], adapower, pendiag = TRUE)
        })

        # 3. Compute common TVP weights (if enabled)
        if (!is.null(initcoefs$common_tvp_effects) && common_tvp_effects) {

          # Note: no pendiag for TVP (time-varying, not autoregressive)
          common_tvp_list <- lapply(seq_along(initcoefs$common_tvp_effects), function(p){
            adaptive_weights(initcoefs$common_tvp_effects[[p]], adapower, pendiag = TRUE)
          })

          # Reformat common TVP weights to match A matrix column structure (block-diagonal by period)
          # For each outcome variable (row): [period1_var1...period1_vard | period2_var1...period2_vard | ...]
          phi_common_weights <- do.call(rbind, lapply(1:n_responses, function(j){
            unlist(lapply(seq_along(common_tvp_list), function(p) common_tvp_list[[p]][j, ]))
          }))

        } else {
          phi_common_weights <- NULL
        }

        # 4. Reformat unique TVP weights (block-diagonal by period)
        # For each subject i, for each period p: weights for [var1, var2, ..., vard]
        # Structure: [period1_vars | period2_vars | ...] for each subject
        phi_unique_weights <- do.call(cbind, lapply(seq_along(tvp_list), function(i){
          # For subject i, stack period blocks: each period has d columns
          do.call(cbind, lapply(seq_along(tvp_list[[i]]), function(p){
            # Period p weights: d x d matrix, use as-is
            tvp_list[[i]][[p]]
          }))
        }))

        # 5. Combine weights in correct order
        if (include_common_effects && !is.null(phi_common_weights)){
          # 4-layer or 3-layer with common base and common TVP
          common_weights <- adaptive_weights(initcoefs$common_effects, adapower, pendiag)
          w_mat <- cbind(common_weights, do.call("cbind", unique_weights_list), phi_common_weights, phi_unique_weights)

        } else if (include_common_effects && is.null(phi_common_weights)){
          # Current default: common base + unique base + unique TVP
          common_weights <- adaptive_weights(initcoefs$common_effects, adapower, pendiag)
          w_mat <- cbind(common_weights, do.call("cbind", unique_weights_list), phi_unique_weights)

        } else if (!include_common_effects && !is.null(phi_common_weights)){
          # No common base, but has common TVP: unique base + common TVP + unique TVP
          w_mat <- cbind(do.call("cbind", unique_weights_list), phi_common_weights, phi_unique_weights)

        } else {
          # No common base, no common TVP: unique base + unique TVP (current common_effects=FALSE)
          w_mat <- cbind(do.call("cbind", unique_weights_list), phi_unique_weights)
        }

      }
      
    }
  }
   
  # Replicate w_mat and scale by ratio parameters
  W <- scale_weights_by_ratios(
    w_mat = w_mat,
    k = k,
    n_pred = n_pred_per_subj[1],
    subgroup = subgroup,
    tvp = tvp,
    include_common_effects = include_common_effects,
    common_tvp_effects = common_tvp_effects,
    ratios_unique = ratios_unique,
    ratios_subgroup = ratios_subgroup,
    ratios_unique_tvp = ratios_unique_tvp,
    ratios_common_tvp = ratios_common_tvp,
    subgroup_membership = subgroup_membership,
    initcoefs = initcoefs
  )

  W
    
}



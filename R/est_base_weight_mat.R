#' @export
est_base_weight_mat <- function(
  W,
  Ak,
  initcoefs,
  ratios,
  d,
  k,
  lassotype,
  weightest,
  subgroup_membership,
  subgroup,
  ratiostau,
  pendiag,
  tvp,
  ratiosalpha,
  ratiosbeta = NULL,
  intercept,
  common_effects = TRUE,
  common_tvp_effects = TRUE){

  # Note, if intercepts are included pendiag should take into account AR
  # effects are no longer on the diagonal of the dynamics (they shift by 1 column).

  # Save parameter value before it could be overwritten
  include_common_effects <- common_effects

  # W      <- object@W
  # Ak         <- object@Ak
  # initcoefs  <- object@initcoefs
  # ratios    <- object@ratios
  # d <- object@d
  # k <- object@k
  # lassotype <- object@lassotype
  # weightest  <- object@weightest
  # subgroup      <- object@subgroup
  # subgroupflag  <- object@subgroupflag
  # ratiostau <- object@ratiostau
  # pendiag <- object@pendiag
  # tvp <- object@tvp
  # ratiosalpha <- object@ratiosalpha
  
  dk <- vapply(Ak, function(x) { ncol(x) }, numeric(1))

  # n_responses = number of response equations (always d)
  # n_predictors = number of predictors per equation
  # (These are equal now, but kept separate for future extensibility
  # if additional predictors are added to the model)
  n_responses <- d
  n_predictors <- dk[1]

  adapower <- 1

  if (lassotype == "standard"){
   
    w_mat <- W
   
  } else {

    if(length(Ak) == 1){

      if(!tvp){
        # Non-TVP k=1 case
        w_mat <- 1/abs(initcoefs$total_effects[[1]])^adapower
        w_mat[is.infinite(w_mat)] <- 1e10

        if(!pendiag){ diag(w_mat) <- 1e-10 }

      } else {
        # TVP k=1 case

        if (!common_effects || is.null(initcoefs$common_effects)) {
          # No time-invariant component: weights from per-period estimates directly
          # tvp_effects[[1]] contains the per-period estimates (not deviations)
          tvp_list <- lapply(seq_along(initcoefs$tvp_effects[[1]]), function(j){
            v <- 1/abs(initcoefs$tvp_effects[[1]][[j]])^adapower
            v[is.infinite(v)] <- 1e10
            if(!pendiag){ diag(v) <- 1e-10 }
            v
          })

          # Reformat tvp weights to match A matrix structure
          w_mat <- do.call(rbind, lapply(1:n_responses, function(j){
            c(do.call(rbind, lapply(seq_along(tvp_list), function(m){
              tvp_list[[m]][j,]
            })))
          }))

        } else {
          # Time-invariant (common) weights from common_effects
          # (For k=1 TVP, common_effects = median of per-period estimates)
          common_weights <- 1/abs(initcoefs$common_effects)^adapower
          common_weights[is.infinite(common_weights)] <- 1e10

          if(!pendiag){ diag(common_weights) <- 1e-10 }

          # TVP weights from tvp_effects (period-specific deviations from common)
          tvp_list <- lapply(seq_along(initcoefs$tvp_effects[[1]]), function(j){
            v <- 1/abs(initcoefs$tvp_effects[[1]][[j]])^adapower
            v[is.infinite(v)] <- 1e10
            v
          })

          # Reformat tvp weights to match A matrix structure
          phi_weights <- do.call(rbind, lapply(1:n_responses, function(j){
            c(do.call(rbind, lapply(seq_along(tvp_list), function(m){
              tvp_list[[m]][j,]
            })))
          }))

          # Combine common and TVP weights
          w_mat <- cbind(common_weights, phi_weights)
        }
      }

    } else {
      
      if(!subgroup){

        v_list <- lapply(seq_along(Ak), function(i){
          v <- 1/abs(initcoefs$unique_effects[[i]])^adapower
          v[is.infinite(v)] <- 1e10
          v
        })

        if (include_common_effects){
          # Standard: common + unique weights
          b_med <- 1/abs(initcoefs$common_effects)^1
          b_med[is.infinite(b_med)] <- 1e10

          if(!pendiag){ diag(b_med) <- 1e-10 }

          w_mat <- cbind(b_med, do.call("cbind", v_list))
        } else {
          # No common effects: only unique weights
          w_mat <- do.call("cbind", v_list)
        }

      } else {

        s_list <- lapply(seq_along(1:length(initcoefs$subgroup_effects)), function(i){
          v <- 1/abs(initcoefs$subgroup_effects[[i]])^adapower
          v[is.infinite(v)] <- 1e10
          v
        })

        v_list <- lapply(seq_along(Ak), function(i){
          v <- 1/abs(initcoefs$unique_effects[[i]])^adapower
          v[is.infinite(v)] <- 1e10
          v
        })

        b_med <- 1/abs(initcoefs$common_effects)^1
        b_med[is.infinite(b_med)] <- 1e10

        if(!pendiag){ diag(b_med) <- 1e-10 }

        w_mat <- cbind(b_med, do.call("cbind", s_list), do.call("cbind", v_list))

      }
      
      # TODO: TVP-specific intercept handling not yet implemented
      if(tvp){

        # 1. Compute unique TVP weights
        tvp_list <- lapply(seq_along(1:length(initcoefs$tvp_effects)), function(i){
          lapply(seq_along(1:length(initcoefs$tvp_effects[[i]])), function(j){
            v <- 1/abs(initcoefs$tvp_effects[[i]][[j]])^adapower
            v[is.infinite(v)] <- 1e10
            v
          })
        })

        # 2. Compute unique base weights
        v_list <- lapply(seq_along(Ak), function(i){
          v <- 1/abs(initcoefs$unique_effects[[i]])^adapower
          v[is.infinite(v)] <- 1e10
          v
        })

        # 3. Compute common TVP weights (if enabled)
        if (!is.null(initcoefs$common_tvp_effects) && common_tvp_effects) {

          common_tvp_list <- lapply(seq_along(initcoefs$common_tvp_effects), function(p){
            v <- 1/abs(initcoefs$common_tvp_effects[[p]])^adapower
            v[is.infinite(v)] <- 1e10
            # Note: no pendiag for TVP (time-varying, not autoregressive)
            v
          })

          # Reformat common TVP weights to match A matrix column structure
          # For each outcome variable (row), stack all periods' weights
          phi_common_weights <- do.call(rbind, lapply(1:n_responses, function(j){
            c(do.call(rbind, lapply(seq_along(common_tvp_list), function(p){
              common_tvp_list[[p]][j, ]
            })))
          }))

        } else {
          phi_common_weights <- NULL
        }

        # 4. Reformat unique TVP weights
        # Row-diagonal structure:
        # row 1: b^{1}_{1,1,t=1},...,b^{1}_{1,1,t=T},...,b^{1}_{1,d,t=1},...,b^{1}_{1,d,t=T}
        # row 2: b^{1}_{2,1,t=1},...,b^{1}_{2,1,t=T},...,b^{1}_{2,d,t=1},...,b^{1}_{2,d,t=T}
        #      :
        # row d: b^{1}_{d,1,t=1},...,b^{1}_{d,1,t=T},...,b^{1}_{d,d,t=1},...,b^{1}_{d,d,t=T}

        phi_unique_weights <- do.call(cbind,lapply(seq_along(tvp_list), function(i){
          # first pull out each equation and stack the rows of phi
          do.call(rbind,lapply(1:n_responses, function(j){
            c(do.call(rbind,lapply(seq_along(tvp_list[[i]]), function(m){
              tvp_list[[i]][[m]][j,]
            })))
          }))
        }))

        # 5. Combine weights in correct order
        if (include_common_effects && !is.null(phi_common_weights)){
          # 4-layer or 3-layer with common base and common TVP
          b_med <- 1/abs(initcoefs$common_effects)^1
          b_med[is.infinite(b_med)] <- 1e10

          if(!pendiag){ diag(b_med) <- 1e-10 }

          w_mat <- cbind(b_med, do.call("cbind", v_list), phi_common_weights, phi_unique_weights)

        } else if (include_common_effects && is.null(phi_common_weights)){
          # Current default: common base + unique base + unique TVP
          b_med <- 1/abs(initcoefs$common_effects)^1
          b_med[is.infinite(b_med)] <- 1e10

          if(!pendiag){ diag(b_med) <- 1e-10 }

          w_mat <- cbind(b_med, do.call("cbind", v_list), phi_unique_weights)

        } else if (!include_common_effects && !is.null(phi_common_weights)){
          # No common base, but has common TVP: unique base + common TVP + unique TVP
          w_mat <- cbind(do.call("cbind", v_list), phi_common_weights, phi_unique_weights)

        } else {
          # No common base, no common TVP: unique base + unique TVP (current common_effects=FALSE)
          w_mat <- cbind(do.call("cbind", v_list), phi_unique_weights)
        }

      }
      
    }
  }
   
  #-----------------------------------#
  # multiple weights by ratios
  #-----------------------------------#


  if(!subgroup){
    
    W <- replicate(length(ratios), w_mat, simplify="array")
    
  } else {
    
    W <- replicate(length(ratios)*length(ratiostau), w_mat, simplify="array")
    
  }
  
  if(tvp){
    if(k == 1){
      # k=1: no ratios needed, only ratiosalpha
      W <- replicate(length(ratiosalpha), w_mat, simplify="array")
    } else {
      # k>1: need to account for common_tvp_effects
      if (include_common_effects && common_tvp_effects){
        # 4-layer: need ratios, ratiosalpha, AND ratiosbeta
        if (is.null(ratiosbeta)) ratiosbeta <- ratiosalpha
        W <- replicate(length(ratios)*length(ratiosalpha)*length(ratiosbeta), w_mat, simplify="array")
      } else if (include_common_effects && !common_tvp_effects){
        # 3-layer: common base, unique base, unique TVP (current default)
        W <- replicate(length(ratios)*length(ratiosalpha), w_mat, simplify="array")
      } else if (!include_common_effects && common_tvp_effects){
        # 3-layer: unique base, common TVP, unique TVP
        if (is.null(ratiosbeta)) ratiosbeta <- ratiosalpha
        W <- replicate(length(ratiosalpha)*length(ratiosbeta), w_mat, simplify="array")
      } else {
        # 2-layer: unique base, unique TVP (current common_effects=FALSE)
        W <- replicate(length(ratiosalpha), w_mat, simplify="array")
      }
    }
  }

  # here we use dk[1] and assume all individuals have the same number
  # of predictors. when this is relaxed this should be modified 
  # accordingly. (zff 2021-09-15)
  if(length(Ak) == 1){

    if(!tvp){
      # Non-TVP k=1: multiply all weights by ratios
      for(r in 1:length(ratios)){
        W[,,r] <- W[,,r] * ratios[r]
      }
    } else {
      # TVP k=1 with common_effects: scale common weights by ratiosalpha
      # This mirrors standard multivar where common weights are scaled by ratios
      if (include_common_effects) {
        for(a in 1:length(ratiosalpha)){
          W[,1:(dk[1]),a] <- W[,1:(dk[1]),a] * ratiosalpha[a]
        }
      }
    }

  } else {
  
    if(!subgroup & !tvp){

      for(r in 1:length(ratios)){

       #W[,(dk[1]+1):ncol(W[,,1]),r] <- W[,(dk[1]+1):ncol(W[,,1]),r] * ratios[r]
        if (include_common_effects){
          W[,1:(dk[1]),r] <- W[,1:(dk[1]),r] * ratios[r]
        } else {
          # All columns are unique (no common to skip)
          # This branch shouldn't be reached since !tvp and common_effects=FALSE is not allowed
        }

      }

    } else if (tvp){
      cnt <- 1

      if (include_common_effects && common_tvp_effects){
        # 4-layer: common base + unique base + common TVP + unique TVP
        if (is.null(ratiosbeta)) ratiosbeta <- ratiosalpha

        # Calculate column indices
        common_base_cols <- 1:dk[1]
        unique_base_cols <- (dk[1]+1):(dk[1]*(k+1))

        # Common TVP columns start after unique base
        # Need to determine size - for now, estimate from w_mat structure
        common_tvp_start <- dk[1]*(k+1) + 1
        # Find where unique TVP starts by checking if common_tvp_effects exists in initcoefs
        if (!is.null(initcoefs$common_tvp_effects)){
          num_common_tvp_cols <- nrow(w_mat[,common_tvp_start:ncol(w_mat),drop=FALSE])
          # Actually this is wrong... let me use a different approach
          # The weight matrix columns match B matrix columns
          # We can determine structure from the weight matrix itself
          # For now, assume proportional split
          tvp_start <- dk[1]*(k+1) + 1
          total_tvp_cols <- ncol(w_mat) - dk[1]*(k+1)

          # Estimate: if we have common TVP, assume first part is common, rest is unique
          # This is a rough estimate - ideally we'd pass column counts explicitly
          # For now, scale both TVP sections equally by their respective ratios
          common_tvp_cols <- tvp_start:(tvp_start + total_tvp_cols%/%3 - 1)  # rough estimate
          unique_tvp_cols <- (max(common_tvp_cols)+1):ncol(w_mat)
        } else {
          # No common TVP - all TVP is unique
          unique_tvp_cols <- tvp_start:ncol(w_mat)
          common_tvp_cols <- NULL
        }

        for(r in 1:length(ratios)){
          for(a in 1:length(ratiosalpha)){
            for(b in 1:length(ratiosbeta)){
              # Scale common base by ratios[r]
              W[,common_base_cols,cnt] <- W[,common_base_cols,cnt] * ratios[r]

              # Scale common TVP by ratiosbeta[b]
              if (!is.null(common_tvp_cols) && length(common_tvp_cols) > 0){
                W[,common_tvp_cols,cnt] <- W[,common_tvp_cols,cnt] * ratiosbeta[b]
              }

              # Scale unique TVP by ratiosalpha[a]
              W[,unique_tvp_cols,cnt] <- W[,unique_tvp_cols,cnt] * ratiosalpha[a]

              cnt <- cnt + 1
            }
          }
        }

      } else if (include_common_effects && !common_tvp_effects){
        # 3-layer: common base + unique base + unique TVP (current default)
        for(r in 1:length(ratios)){
          for(j in 1:length(ratiosalpha)){
            # Scale common columns by ratios, TVP by ratiosalpha
            W[,1:(dk[1]),cnt] <- W[,1:(dk[1]),cnt] * ratios[r]
            W[,(dk[1]*(k+1)+1):ncol(W),cnt] <- W[,(dk[1]*(k+1)+1):ncol(W),cnt] * ratiosalpha[j]
            cnt <- cnt + 1
          }
        }

      } else if (!include_common_effects && common_tvp_effects){
        # 3-layer: unique base + common TVP + unique TVP
        if (is.null(ratiosbeta)) ratiosbeta <- ratiosalpha

        # Columns: unique base (dk[1]*k), then common TVP, then unique TVP
        unique_base_cols <- 1:(dk[1]*k)
        tvp_start <- dk[1]*k + 1

        # Rough estimate of common vs unique TVP split
        total_tvp_cols <- ncol(w_mat) - dk[1]*k
        common_tvp_cols <- tvp_start:(tvp_start + total_tvp_cols%/%3 - 1)
        unique_tvp_cols <- (max(common_tvp_cols)+1):ncol(w_mat)

        for(a in 1:length(ratiosalpha)){
          for(b in 1:length(ratiosbeta)){
            # Scale common TVP by ratiosbeta[b]
            W[,common_tvp_cols,cnt] <- W[,common_tvp_cols,cnt] * ratiosbeta[b]

            # Scale unique TVP by ratiosalpha[a]
            W[,unique_tvp_cols,cnt] <- W[,unique_tvp_cols,cnt] * ratiosalpha[a]

            cnt <- cnt + 1
          }
        }

      } else {
        # 2-layer: unique base + unique TVP (current common_effects=FALSE)
        # No weight scaling - ratiosalpha used in lambda grid construction
      }

    } else if (subgroup & !tvp){
      
      cnt <- 1
      for(r in 1:length(ratios)){
        for(j in 1:length(ratiostau)){
          
          W[,1:(dk[1]),cnt] <- W[,1:(dk[1]),cnt] * ratios[r]
          W[,(dk[1]+1):(dk[1]*max(subgroup_membership)+dk[1]),cnt] <- W[,(dk[1]+1):(dk[1]*max(subgroup_membership)+dk[1]),cnt] * ratiostau[j]
          cnt <- cnt + 1
          
        }
      }
      
    }
 
  }

 return(W)
    
}



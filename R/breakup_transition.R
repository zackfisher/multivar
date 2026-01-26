breakup_transition <- function(B, Ak, ndk, intercept, thresh, subgroup_membership, subgroup, tvp, ntk, breaks, common_effects = TRUE, common_tvp_effects = TRUE){

  # Save parameter value before potential variable reuse
  include_common_effects <- common_effects

  # B <- fit[[1]][,,which.min(colMeans(fit[[2]]))]
  # Ak <- object@Ak
  # ndk <- object@ndk
  # intercept <- object@intercept
  # thresh <- object@thresh
  # subgroup <- object@subgroup
  # subgroupflag <- object@subgroupflag
  # tvp <- object@tvp
  # ntk <- object@ntk
  # breaks <- object@breaks

  #if(!intercept){ B <- B[,-1] }
  if(intercept){
   ndk <- vapply(Ak, function(x) { ncol(x) }, numeric(1))
  }
  
  if(length(Ak) == 1){

    if(!tvp){
      # Non-TVP k=1 case
      # Set column and row names properly (especially for intercept handling)
      colnames(B) <- colnames(Ak[[1]])
      if(intercept){
        rownames(B) <- colnames(Ak[[1]])[-1]  # exclude Intercept from rownames
      } else {
        rownames(B) <- colnames(Ak[[1]])
      }

      common_mat  <- B
      unique_mats <- list(B)
      total_mats  <- list(B)
      diff_mats   <- NULL
      subgrp_mats <- NULL
      tvp_mats    <- NULL

    } else {
      # TVP k=1 case
      # Base effects (columns 1:ndk[1])
      base_mat <- B[, 1:ndk[1], drop=FALSE]
      colnames(base_mat) <- colnames(Ak[[1]])
      if(intercept){
        rownames(base_mat) <- colnames(Ak[[1]])[-1]
      } else {
        rownames(base_mat) <- colnames(Ak[[1]])
      }

      common_mat  <- base_mat
      unique_mats <- list(base_mat)

      # TVP effects (columns (ndk[1]+1):ncol(B))
      # Extract TVP coefficients for each period
      tvp_start <- ndk[1] + 1
      num_periods <- length(breaks[[1]])

      tvp_mats <- list(
        lapply(1:num_periods, function(p){
          # For period p, extract columns for all variables
          period_cols <- seq(tvp_start + (p-1), by=num_periods, length.out=ndk[1])
          mat <- B[, period_cols, drop=FALSE]
          colnames(mat) <- colnames(Ak[[1]])
          if(intercept){
            rownames(mat) <- colnames(Ak[[1]])[-1]
          } else {
            rownames(mat) <- colnames(Ak[[1]])
          }
          mat
        })
      )

      # Total effects for each period = base + tvp
      total_mats <- list(
        lapply(1:num_periods, function(p){
          g <- base_mat + tvp_mats[[1]][[p]]
          g[abs(g) < thresh] <- 0
          g
        })
      )

      diff_mats   <- NULL
      subgrp_mats <- NULL
    }

  } else {

    # common mat (only extract if common_effects = TRUE)
    if (include_common_effects){
      # indices assume intercept has been removed
      first_com_col_index   <- 1 # ifelse(intercept, 1, 2)
      final_com_col_index   <- ndk[1]

      common_mat <- B[,first_com_col_index:final_com_col_index]
      # colnames = predictors (includes Intercept if present)
      # rownames = outcomes (d variables, no intercept)
      colnames(common_mat) <- colnames(Ak[[1]])
      if(intercept){
        rownames(common_mat) <- colnames(Ak[[1]])[-1]  # exclude Intercept from rownames
      } else {
        rownames(common_mat) <- colnames(Ak[[1]])
      }
    } else {
      # No common effects
      common_mat <- NULL
    }
    
    if(!subgroup & !tvp){

      if (include_common_effects){
        # Standard: common + unique columns
        first_ind_col_indices <- cumsum(ndk) + 1
        final_ind_col_indices <- first_ind_col_indices + ndk - 1

        # unique mats
        unique_mats <- lapply(seq_along(ndk), function(i){
          mat <- B[,first_ind_col_indices[i]:final_ind_col_indices[i]]
          colnames(mat) <- colnames(Ak[[i]])
          if(intercept){
            rownames(mat) <- colnames(Ak[[i]])[-1]
          } else {
            rownames(mat) <- colnames(Ak[[i]])
          }
          mat
        })

        # total mats
        total_mats <- lapply(unique_mats, function(mat){
          g <- mat + common_mat;
          g[abs(g) < thresh] <- 0
          g
        })

      } else {
        # No common effects: only unique columns
        # Column indices shift (no common block to skip)
        first_ind_col_indices <- c(1, cumsum(ndk)[-length(ndk)] + 1)
        final_ind_col_indices <- cumsum(ndk)

        # unique mats
        unique_mats <- lapply(seq_along(ndk), function(i){
          mat <- B[,first_ind_col_indices[i]:final_ind_col_indices[i]]
          colnames(mat) <- colnames(Ak[[i]])
          if(intercept){
            rownames(mat) <- colnames(Ak[[i]])[-1]
          } else {
            rownames(mat) <- colnames(Ak[[i]])
          }
          mat
        })

        # total mats = unique mats (no common to add)
        total_mats <- lapply(unique_mats, function(mat){
          g <- mat
          g[abs(g) < thresh] <- 0
          g
        })
      }

      subgrp_mats <- NULL
      tvp_mats    <- NULL
      
    } else if (subgroup & !tvp) {
      
      first_sub_col_indices <- cumsum(rep(ndk[1],max(subgroup_membership))) + 1
      final_sub_col_indices <- first_sub_col_indices + ndk[1] - 1
      
      first_ind_col_indices <- cumsum(ndk) + (ndk[1]*max(subgroup_membership)) + 1 
      final_ind_col_indices <- first_ind_col_indices + ndk - 1
      
      # subgrp mats
      subgrp_mats <- lapply(seq_along(ndk), function(i){
        mat <- B[,first_sub_col_indices[subgroup_membership[i]]:final_sub_col_indices[subgroup_membership[i]]]
        colnames(mat) <- colnames(Ak[[i]])
        if(intercept){
          rownames(mat) <- colnames(Ak[[i]])[-1]
        } else {
          rownames(mat) <- colnames(Ak[[i]])
        }
        mat
      })

      # unique mats
      unique_mats <- lapply(seq_along(ndk), function(i){
        mat <- B[,first_ind_col_indices[i]:final_ind_col_indices[i]]
        colnames(mat) <- colnames(Ak[[i]])
        if(intercept){
          rownames(mat) <- colnames(Ak[[i]])[-1]
        } else {
          rownames(mat) <- colnames(Ak[[i]])
        }
        mat
      })
      
      # total mats
      total_mats <- lapply(seq_along(unique_mats), function(i){
        g <- unique_mats[[i]] + subgrp_mats[[i]] + common_mat;
        g[abs(g) < thresh] <- 0
        g
      })
      
      tvp_mats    <- NULL
      
    } else if (tvp){

      if (include_common_effects){
        # Standard: common + unique + tvp columns
        first_ind_col_indices <- cumsum(ndk) + 1
        final_ind_col_indices <- first_ind_col_indices + ndk - 1
      } else {
        # No common: unique + tvp columns (shift indices)
        first_ind_col_indices <- c(1, cumsum(ndk)[-length(ndk)] + 1)
        final_ind_col_indices <- cumsum(ndk)
      }

      # unique mats
      unique_mats <- lapply(seq_along(ndk), function(i){
        mat <- B[,first_ind_col_indices[i]:final_ind_col_indices[i]]
        colnames(mat) <- colnames(Ak[[i]])
        if(intercept){
          rownames(mat) <- colnames(Ak[[i]])[-1]
        } else {
          rownames(mat) <- colnames(Ak[[i]])
        }
        mat
      })

      # Extract common TVP effects (if enabled)
      if (common_tvp_effects) {

        num_periods <- length(breaks[[1]])
        common_tvp_start <- max(final_ind_col_indices) + 1
        common_tvp_ncols <- num_periods * ndk[1]
        common_tvp_end <- common_tvp_start + common_tvp_ncols - 1

        common_tvp_B <- B[, common_tvp_start:common_tvp_end, drop=FALSE]

        # Extract per-period common TVP matrices
        common_tvp_mats <- lapply(1:num_periods, function(p){
          period_cols <- seq((p-1)*ndk[1] + 1, p*ndk[1])
          mat <- common_tvp_B[, period_cols, drop=FALSE]
          colnames(mat) <- colnames(Ak[[1]])
          if(intercept){
            rownames(mat) <- colnames(Ak[[1]])[-1]
          } else {
            rownames(mat) <- colnames(Ak[[1]])
          }
          mat
        })

        unique_tvp_start <- common_tvp_end + 1

      } else {
        common_tvp_mats <- NULL
        unique_tvp_start <- max(final_ind_col_indices) + 1
      }

      # unique TVP mats
      # row 1: b^{1}_{1,1,t=1},...,b^{1}_{1,1,t=T},...,b^{1}_{1,d,t=1},...,b^{1}_{1,d,t=T}
      # row 2: b^{1}_{2,1,t=1},...,b^{1}_{2,1,t=T},...,b^{1}_{2,d,t=1},...,b^{1}_{2,d,t=T}
      #      :
      # row d: b^{1}_{d,1,t=1},...,b^{1}_{d,1,t=T},...,b^{1}_{d,d,t=1},...,b^{1}_{d,d,t=T}

      cols_k <- unlist(lapply(breaks,function(g){length(g)*ndk[1]}))

      first_tvp_col_indices <- c(0, cumsum(cols_k)[-length(Ak)]) + unique_tvp_start
      final_tvp_col_indices <- cumsum(cols_k) + unique_tvp_start - 1

      # Unique TVP extraction
      unique_tvp_mats <- lapply(seq_along(unique_mats), function(i){

        mat <- B[,first_tvp_col_indices[i]:final_tvp_col_indices[i]]

        tvp_mats_inner <- lapply(1:nrow(mat), function(j){
          m <- matrix(c(mat[j,]), nrow = ndk[1], byrow = TRUE)
          g <- unlist(lapply(breaks[[i]],function(g){length(g)}))
          m[,rep(1:ncol(m), times = g)]
        })

        lapply(1:ntk[1], function(j){
          eq_d <- lapply(1:ndk[1], function(eq){
            tvp_mats_inner[[eq]][,j]
          })
          do.call(rbind, eq_d)
        })

      })

      # Total effects: combine all layers
      total_mats <- lapply(seq_along(unique_tvp_mats), function(i){
        lapply(seq_along(unique_tvp_mats[[i]]), function(t){

          # Determine which period time point t is in
          period_idx <- which(sapply(breaks[[i]], function(b) t %in% b))

          # Build total: base + TVP
          if (include_common_effects){
            total <- common_mat + unique_mats[[i]]
          } else {
            total <- unique_mats[[i]]
          }

          # Add common TVP (if enabled)
          if (common_tvp_effects && !is.null(common_tvp_mats)){
            total <- total + common_tvp_mats[[period_idx]]
          }

          # Add unique TVP
          total <- total + unique_tvp_mats[[i]][[t]]

          # Threshold
          total[abs(total) < thresh] <- 0
          total
        })
      })

      subgrp_mats <- NULL
      
    }
    
  }


  # For backward compatibility, assign unique_tvp_mats to tvp_mats
  # (tvp in the return list refers to unique TVP effects)
  if (exists("unique_tvp_mats")) {
    tvp_mats <- unique_tvp_mats
  }

  res <- list(
    common = common_mat,
    subgrp = subgrp_mats,
    unique = unique_mats,
    tvp  = tvp_mats,
    common_tvp = if(exists("common_tvp_mats")) common_tvp_mats else NULL,
    total  = total_mats
  )
  return(res)
}
  
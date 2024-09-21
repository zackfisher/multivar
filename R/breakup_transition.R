breakup_transition <- function(B, Ak, ndk, intercept, thresh, subgroup_membership, subgroup, tvp, ntk, breaks){
  
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
  
  if(!intercept){ B <- B[,-1] } 
  
  if(length(Ak) == 1){
    
    common_mat  <- B
    unique_mats <- list(B)
    total_mats  <- list(B)
    diff_mats   <- NULL
    
  } else {
    
    # indices assume intercept has been removed
    first_com_col_index   <- 1 # ifelse(intercept, 1, 2)
    final_com_col_index   <- ndk[1] 
    
    # common mat
    common_mat <- B[,first_com_col_index:final_com_col_index]
    rownames(common_mat) <- colnames(common_mat) <- colnames(Ak[[1]])
    
    if(!subgroup & !tvp){
      
      first_ind_col_indices <- cumsum(ndk) + 1
      final_ind_col_indices <- first_ind_col_indices + ndk - 1
      
      # unique mats
      unique_mats <- lapply(seq_along(ndk), function(i){
        mat <- B[,first_ind_col_indices[i]:final_ind_col_indices[i]]
        rownames(mat) <- colnames(mat) <- colnames(Ak[[i]])
        mat
      })
      
      
      # total mats
      total_mats <- lapply(unique_mats, function(mat){
        g <- mat + common_mat;
        g[abs(g) < thresh] <- 0
        g
      })
      
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
        rownames(mat) <- colnames(mat) <- colnames(Ak[[i]])
        mat
      })
      
      # unique mats
      unique_mats <- lapply(seq_along(ndk), function(i){
        mat <- B[,first_ind_col_indices[i]:final_ind_col_indices[i]]
        rownames(mat) <- colnames(mat) <- colnames(Ak[[i]])
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
      
      first_ind_col_indices <- cumsum(ndk) + 1
      final_ind_col_indices <- first_ind_col_indices + ndk - 1
      
      # unique mats
      unique_mats <- lapply(seq_along(ndk), function(i){
        mat <- B[,first_ind_col_indices[i]:final_ind_col_indices[i]]
        rownames(mat) <- colnames(mat) <- colnames(Ak[[i]])
        mat
      })
      
      # tvp mats
              # row 1: b^{1}_{1,1,t=1},...,b^{1}_{1,1,t=T},...,b^{1}_{1,d,t=1},...,b^{1}_{1,d,t=T}
        # row 2: b^{1}_{2,1,t=1},...,b^{1}_{2,1,t=T},...,b^{1}_{2,d,t=1},...,b^{1}_{2,d,t=T}
        #      :
        # row d: b^{1}_{d,1,t=1},...,b^{1}_{d,1,t=T},...,b^{1}_{d,d,t=1},...,b^{1}_{d,d,t=T}
      
      #splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))
      #bandwidth_breaks <- seq(from = bandwidth, to = ntk[1], by = bandwidth)
      #break_list <- splitAt(c(1:ntk[1]), bandwidth_breaks)
      
      #first_tvp_col_indices <- c(0, cumsum(rep(sum(object@ntk),ndk[1]))[-ndk[1]]) + max(final_ind_col_indices) + 1
      #final_tvp_col_indices <- cumsum(rep(sum(object@ntk),d)) + max(final_ind_col_indices) 
      cols_k <- unlist(lapply(breaks,function(g){length(g)*ndk[1]}))
      #cols_k_old <- rep(sum(length(break_list)*ndk[1]),length(ntk))
      
      first_tvp_col_indices <- c(0, cumsum(cols_k)[-length(Ak)]) + max(final_ind_col_indices) + 1
      final_tvp_col_indices <- cumsum(cols_k) + max(final_ind_col_indices) 
      
      # total mats
      tvp_mats <- lapply(seq_along(unique_mats), function(i){
        
        mat <- B[,first_tvp_col_indices[i]:final_tvp_col_indices[i]]
        
        tvp_mats <- lapply(1:nrow(mat), function(j){
          m <- matrix(c(mat[j,]), nrow = ndk[1], byrow = TRUE)
          g <- unlist(lapply(breaks[[i]],function(g){length(g)}))
          m[,rep(1:ncol(m), times = g)]
        })
        
        lapply(1:ntk[1], function(j){
          eq_d <- lapply(1:ndk[1], function(eq){
            tvp_mats[[eq]][,j]
          })
          do.call(rbind, eq_d)
        })
        
      })

      # total mats
      total_mats <- lapply(seq_along(tvp_mats), function(i){
        lapply(seq_along(tvp_mats[[i]]), function(j){
          g <- common_mat + unique_mats[[i]] + tvp_mats[[i]][[j]]
          g[abs(g) < thresh] <- 0
          g
        })
      })
      
      subgrp_mats <- NULL
      
    }
    
  }
  

  res <- list(
    common = common_mat,
    subgrp = subgrp_mats,
    unique = unique_mats,
    tvp  = tvp_mats,
    total  = total_mats
  )
  return(res)
}
  
breakup_transition <- function(B, Ak, ndk, intercept, thresh, subgroup,subgroupflag){
  
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
    
    if(!subgroupflag){
      
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
      
    } else {
      
      first_sub_col_indices <- cumsum(rep(ndk[1],max(subgroup))) + 1
      final_sub_col_indices <- first_sub_col_indices + ndk[1] - 1
      
      first_ind_col_indices <- cumsum(ndk) + (ndk[1]*max(subgroup)) + 1 
      final_ind_col_indices <- first_ind_col_indices + ndk - 1
      
      # subgrp mats
      subgrp_mats <- lapply(seq_along(ndk), function(i){
        mat <- B[,first_sub_col_indices[subgroup[i]]:final_sub_col_indices[subgroup[i]]]
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
    }
    
  }
  

  res <- list(
    common = common_mat,
    subgrp = subgrp_mats,
    unique = unique_mats,
    total  = total_mats
  )
  return(res)
}
  
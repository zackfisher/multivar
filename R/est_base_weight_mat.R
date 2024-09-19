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
  subgroup,     
  subgroupflag, 
  ratiostau,
  pendiag,
  tvp,
  ratiosalpha){
  
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
  
  adapower <- 1
  
  if (lassotype == "standard"){
   
    w_mat <- W
   
  } else {

    if(length(Ak) == 1){
      
      w_mat <- 1/abs(initcoefs$total_effects[[1]])^adapower
      w_mat[is.infinite(w_mat)] <- 1e10
      if(!pendiag){ diag(w_mat) <- 1e-10 }
      
    } else {
      
      if(!subgroupflag){
        
        v_list <- lapply(seq_along(Ak), function(i){
          v <- 1/abs(initcoefs$unique_effects[[i]])^adapower
          v[is.infinite(v)] <- 1e10
          v
        })
      
        b_med <- 1/abs(initcoefs$common_effects)^1
        b_med[is.infinite(b_med)] <- 1e10
        if(!pendiag){ diag(b_med) <- 1e-10 }
      
        w_mat <- cbind(b_med, do.call("cbind", v_list))
        
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
        
        w_mat <- cbind(b_med, do.call("cbind", s_list),do.call("cbind", v_list))

      }
      
      if(tvp){
        
        tvp_list <- lapply(seq_along(1:length(initcoefs$tvp_effects)), function(i){
          lapply(seq_along(1:length(initcoefs$tvp_effects[[i]])), function(j){
            v <- 1/abs(initcoefs$tvp_effects[[i]][[j]])^adapower
            v[is.infinite(v)] <- 1e10
            v
          })
        })
        
        v_list <- lapply(seq_along(Ak), function(i){
          v <- 1/abs(initcoefs$unique_effects[[i]])^adapower
          v[is.infinite(v)] <- 1e10
          v
        })
        
        b_med <- 1/abs(initcoefs$common_effects)^1
        b_med[is.infinite(b_med)] <- 1e10
        if(!pendiag){ diag(b_med) <- 1e-10 }
        
        # reformat tvp_list
        # it should look like a row-diagonal matrix
        #
        # y1 (row 1): k1,d1, k1,d2, ..., k2,d1.
        
        # row 1: b^{1}_{1,1,t=1},...,b^{1}_{1,1,t=T},...,b^{1}_{1,d,t=1},...,b^{1}_{1,d,t=T}
        # row 2: b^{1}_{2,1,t=1},...,b^{1}_{2,1,t=T},...,b^{1}_{2,d,t=1},...,b^{1}_{2,d,t=T}
        #      :
        # row d: b^{1}_{d,1,t=1},...,b^{1}_{d,1,t=T},...,b^{1}_{d,d,t=1},...,b^{1}_{d,d,t=T}
        
        phi_weights <- do.call(cbind,lapply(seq_along(tvp_list), function(i){
          # first pull out each equation and stack the rows of phi 
          do.call(rbind,lapply(c(1:d[1]), function(j){
            c(do.call(rbind,lapply(seq_along(tvp_list[[i]]), function(m){
              tvp_list[[i]][[m]][j,]
            })))
          }))
        }))

        w_mat <- cbind(b_med, do.call("cbind", v_list), phi_weights)
        
      }
      
    }
  }
   
  #-----------------------------------#
  # multiple weights by ratios
  #-----------------------------------#


  if(!subgroupflag){
    
    W <- replicate(length(ratios), w_mat, simplify="array")
    
  } else {
    
    W <- replicate(length(ratios)*length(ratiostau), w_mat, simplify="array")
    
  }
  
  if(tvp){
    W <- replicate(length(ratios)*length(ratiosalpha), w_mat, simplify="array")
  }

  # here we use d[1] and assume all individuals have the same number
  # of predictors. when this is relaxed this should be modified 
  # accordingly. (zff 2021-09-15)
  if(length(Ak) == 1){
    
    for(r in 1:length(ratios)){
      
      W[,,r] <- W[,,r] * ratios[r]
      
    }
    
  } else {
  
    if(!subgroupflag & !tvp){
      
      for(r in 1:length(ratios)){
        
       #W[,(d[1]+1):ncol(W[,,1]),r] <- W[,(d[1]+1):ncol(W[,,1]),r] * ratios[r]
        W[,1:(d[1]),r] <- W[,1:(d[1]),r] * ratios[r]
        
      }
      
    } else if (tvp){
      cnt <- 1
      for(r in 1:length(ratios)){
        for(j in 1:length(ratiosalpha)){
          
          W[,1:(d[1]),cnt] <- W[,1:(d[1]),cnt] * ratios[r]
          W[,(d[1]*(k+1)+1):ncol(W),cnt] <- W[,(d[1]*(k+1)+1):ncol(W),cnt] * ratiosalpha[j]
          cnt <- cnt + 1
          
        }
      }
      
    } else if (subgroupflag & !tvp){
      
      cnt <- 1
      for(r in 1:length(ratios)){
        for(j in 1:length(ratiostau)){
          
          W[,1:(d[1]),cnt] <- W[,1:(d[1]),cnt] * ratios[r]
          W[,(d[1]+1):(d[1]*max(subgroup)+d[1]),cnt] <- W[,(d[1]+1):(d[1]*max(subgroup)+d[1]),cnt] * ratiostau[j]
          cnt <- cnt + 1
          
        }
      }
      
    }
 
  }

 return(W)
    
}



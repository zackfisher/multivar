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
  ratiostau){
  
  adapower <- 1
  
  if (lassotype == "standard"){
   
    w_mat <- W
   
  } else {

    if(length(Ak) == 1){
      
      w_mat <- 1/abs(initcoefs$total_effects[[1]])^adapower
      w_mat[is.infinite(w_mat)] <- 1e10
      
    } else {
      
      if(!subgroupflag){
        
        v_list <- lapply(seq_along(Ak), function(i){
          v <- 1/abs(initcoefs$unique_effects[[i]])^adapower
          v[is.infinite(v)] <- 1e10
          v
        })
      
        b_med <- 1/abs(initcoefs$common_effects)^1
        b_med[is.infinite(b_med)] <- 1e10
      
        w_mat <- cbind(b_med, do.call("cbind", v_list))
        
      } else {
        
        
        s_list <- lapply(seq_along(1:length(initcoefs$subgroup_effects)), function(i){
          v <- 1/abs(initcoefs$subgroup_effects[[i]])^adapower
          v[is.infinite(v)] <- 1e10
          v
        })
        
        v_list <- lapply(seq_along(Ak), function(i){
          v <- 1/abs(initcoefs$unique_effects)^adapower
          v[is.infinite(v)] <- 1e10
          v
        })
        
        b_med <- 1/abs(initcoefs$common_effects)^1
        b_med[is.infinite(b_med)] <- 1e10
        w_mat <- cbind(b_med, do.call("cbind", s_list),do.call("cbind", v_list))

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

  # here we use d[1] and assume all individuals have the same number
  # of predictors. when this is relaxed this should be modified 
  # accordingly. (zff 2021-09-15)
  if(length(Ak) == 1){
    
    for(r in 1:length(ratios)){
      
      W[,,r] <- W[,,r] * ratios[r]
      
    }
    
  } else {
  
    if(!subgroupflag){
      
      for(r in 1:length(ratios)){
        
       #W[,(d[1]+1):ncol(W[,,1]),r] <- W[,(d[1]+1):ncol(W[,,1]),r] * ratios[r]
        W[,1:(d[1]),r] <- W[,1:(d[1]),r] * ratios[r]
        
      }
      
    } else {
      
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



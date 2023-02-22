#' @export
est_base_weight_mat <- function(
  W,
  Ak, 
  bk,
  ratios,
  d, 
  k, 
  lassotype, 
  weightest,
  subgroup,
  ratiostau){
  
  # W<-object@W
  # Ak<-object@Ak
  # bk<- object@bk
  # ratios<-object@ratios 
  # d<-object@d 
  # k<-object@k 
  # lassotype<-object@lassotype 
  # weightest<-object@weightest
  # subgroup <- -object@subgroup
  
  # eventually make this an argument
  adapower <- 1
  
 if (lassotype == "standard"){
   
   # w_mat <- matrix(1, nrow = ncol(dat@bk[[1]]), ncol = ncol(dat@A))
   w_mat <- W
   
 } else {
   
   #-----------------------------------#
   # estimate k= 1,..,K transition mats
   #-----------------------------------#
   
   if (weightest == "ols") {
     
     b_w <- lapply(seq_along(Ak), function(z){
        b_ols_i <- matrix(0, ncol(Ak[[z]]),ncol(Ak[[z]]))
         for(i in 1:ncol(Ak[[z]])){
          for(j in 1:ncol(Ak[[z]])){
           b_ols_i[j,i] <- lm.fit(as.matrix(Ak[[z]])[,i,drop=F],as.matrix(bk[[z]])[,j,drop=F])$coefficients
          }
         }
       b_ols_i
     })
     
   } else if (weightest == "lasso") {
     
     b_w  <- lapply(seq_along(Ak),function(g){
        lasso.fit <- glmnet::cv.glmnet(
          diag(ncol(Ak[[g]])) %x% Ak[[g]], 
          as.vector(bk[[g]]), 
          family = "gaussian", 
          alpha = 1, 
          standardize = FALSE, 
          nfolds = 10
        )
        t(matrix(as.vector(predict(
          lasso.fit, 
          diag(ncol(Ak[[g]])) %x% Ak[[g]], 
          type="coefficients", 
          s="lambda.1se"
        ))[-1], ncol(Ak[[g]]), ncol(Ak[[g]])))
     })
     
   } else if (weightest == "ridge") {
     
     b_w  <- lapply(seq_along(Ak),function(g){
        lasso.fit <- glmnet::cv.glmnet(
          diag(ncol(Ak[[g]])) %x% Ak[[g]], 
          as.vector(bk[[g]]), 
          family = "gaussian", 
          alpha = 0, 
          standardize = FALSE, 
          nfolds = 10
        )
        t(matrix(as.vector(predict(
          lasso.fit, 
          diag(ncol(Ak[[g]])) %x% Ak[[g]], 
          type="coefficients", 
          s="lambda.1se"
        ))[-1], ncol(Ak[[g]]), ncol(Ak[[g]])))
     })
     
   } else if (weightest == "var") {
     
     b_w  <- lapply(seq_along(Ak),function(g){
        fit<- vars::VAR(as.matrix(Ak[[g]]), p=1, type="none")$varresult
        as.matrix(do.call("rbind",lapply(seq_along(colnames(Ak[[g]])), function(x) {
           fit[[x]]$coefficients
        })))
     })
     
   } else if (weightest == "variable") {
     
      if(ncol(Ak[[1]]) >= nrow(Ak[[1]])){
    
        b_w <- lapply(seq_along(Ak), function(z){
          b_ols_i <- matrix(0, ncol(Ak[[z]]),ncol(Ak[[z]]))
           for(i in 1:ncol(Ak[[z]])){
            for(j in 1:ncol(Ak[[z]])){
             b_ols_i[j,i] <- lm.fit(as.matrix(Ak[[z]])[,i,drop=F],as.matrix(bk[[z]])[,j,drop=F])$coefficients
            }
           }
         b_ols_i
       })
        
      } else {
    
        b_w  <- lapply(seq_along(Ak),function(g){
           fit<- vars::VAR(as.matrix(Ak[[g]]), p=1, type="none")$varresult
           as.matrix(do.call("rbind",lapply(seq_along(colnames(Ak[[g]])), function(x) {
              fit[[x]]$coefficients
           })))
        })
    
       }

   }
 
   
   
   
   #-----------------------------------#
   # 
   #-----------------------------------#
   
   # convert list of coef matrices (b_w) to array
   a <- array(unlist(b_w), c(dim(b_w[[1]]), length(b_w)))
    
    if(length(Ak) == 1){
      
      v_list <- lapply(seq_along(Ak), function(i){
        v <- 1/abs(b_w[[i]])^adapower
        v[is.infinite(v)] <- 1e10
        v
      })
        
      w_mat <- v_list[[1]]
      
      
    } else {
      
      # create matrix of element-wise medians for group effects
      b_med <- apply(a, 1:2, median)
      
      # create array of element-wise medians for subgroup effects
      if(is.null(subgroup)){
        
        v_list <- lapply(seq_along(Ak), function(i){
          v <- 1/abs(b_w[[i]] - b_med)^adapower
          v[is.infinite(v)] <- 1e10
          v
        })

        b_med <- 1/abs(b_med)^1
        b_med[is.infinite(b_med)] <- 1e10
        
        w_mat <- cbind(b_med, do.call("cbind", v_list))
        
      } else {
        
        # subgroup level weights
        b_med_subgrp <- lapply(seq_along(Ak), function(i){
          apply(a[,,which(subgroup==subgroup[i])], 1:2, median)
        })
        
        s_list <- lapply(seq_along(Ak), function(i){
          v <- 1/abs(b_med_subgrp[[i]] - b_med)^adapower
          v[is.infinite(v)] <- 1e10
          v
        })
        
        # individual level weights
        v_list <- lapply(seq_along(Ak), function(i){
          v <- 1/abs(b_w[[i]] - b_med - b_med_subgrp[[i]])^adapower
          v[is.infinite(v)] <- 1e10
          v
        })
        
        b_med <- 1/abs(b_med)^1
        b_med[is.infinite(b_med)] <- 1e10
        w_mat <- cbind(b_med, do.call("cbind", s_list),do.call("cbind", v_list))
        
        
        
      }
      

      
    }

 }
  
  if(is.null(subgroup)){
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
    
    if(is.null(subgroup)){
      for(r in 1:length(ratios)){
       #W[,(d[1]+1):ncol(W[,,1]),r] <- W[,(d[1]+1):ncol(W[,,1]),r] * ratios[r]
        W[,1:(d[1]),r] <- W[,1:(d[1]),r] * ratios[r]
      }
    } else {
      for(r in 1:length(ratios)){
        for(j in 1:length(ratios)){
        W[,1:(d[1]),r] <- W[,1:(d[1]),r] * ratios[r]
        W[,(d[1]+1):(d[1]*max(subgroup)+d[1]),j] <- W[,(d[1]+1):(d[1]*max(subgroup)+d[1]),j] * ratiostau[r]
        }
      }
    }
    
  }

 return(W)
    
}

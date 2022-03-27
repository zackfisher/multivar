#' @export
est_base_weight_mat <- function(
  W,
  Ak, 
  bk,
  ratios,
  d, 
  k, 
  lassotype, 
  weightest){
  
  # W <- object@W
  # Ak <- object@Ak
  # bk <-   object@bk
  # ratios <-  object@ratios 
  # d <-   object@d
  # k <-  object@k
  # lassotype <-  object@lassotype 
  # weightest <-  object@weightest
  
 if (lassotype == "standard"){
   
   # w_mat <- matrix(1, nrow = ncol(dat@bk[[1]]), ncol = ncol(dat@A))
   w_mat <- W
   
 } else {
   
   # if(weightest == "ridge"){
   #   
   #   b_w  <- lapply(seq_along(Ak),function(g){
   #     A_multivar_big <- as.matrix(Ak[[g]])
   #     b_multivar_big <- as.matrix(bk[[g]])
   #     A   <- as.matrix(Matrix::bdiag(replicate(d, A_multivar_big, simplify=FALSE)))
   #     b   <- c(b_multivar_big)
   #     fit <- glmnet::cv.glmnet(A,b,intercept=FALSE,standardize=FALSE,alpha=0)
   #     b_ols_i <- t(matrix(c(as.matrix(coef(fit,s = "lambda.min")))[-1], d, d, byrow=F))
   #     b_ols_i
   #   })
   #   
   # } else if (weightest == "ols") {
   
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
     
   } else if (weightest == "VAR") {
     
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
 
   
    a <- array(unlist(b_w), c(dim(b_w[[1]]), length(b_w)))
    
    if(length(Ak) == 1){
      
      v_list <- lapply(seq_along(Ak), function(i){
        v <- 1/abs(b_w[[i]])^1
        v[is.infinite(v)] <- 1e10
        v
      })
        
      w_mat <- v_list[[1]]
      
      
    } else {
      
      b_med <- apply(a, 1:2, median)
      
      v_list <- lapply(seq_along(Ak), function(i){
        v <- 1/abs(b_w[[i]] - b_med)^1
        v[is.infinite(v)] <- 1e10
        v
      })
      
      b_med <- 1/abs(b_med)^1
      b_med[is.infinite(b_med)] <- 1e10
        
      w_mat <- cbind(b_med, do.call("cbind", v_list))
      
      if(pendiff){
        
        d_list <- lapply(seq_along(Ak), function(i){
          d <- 1/abs(b_w[[i]])^1
          d[is.infinite(d)] <- 1e10
          d
        })
        w_mat <- cbind(w_mat, do.call("cbind", d_list))
        
      }
      
    }

  }
  
  W <- replicate(length(ratios), w_mat, simplify="array")
  
  # here we use d[1] and assume all individuals have the same number
  # of predictors. when this is relaxed this should be modified 
  # accordingly. (zff 2021-09-15)
  if(length(Ak) == 1){
    
    for(r in 1:length(ratios)){
      W[,,r] <- W[,,r] * ratios[r]
    }
    
  } else {
    
    for(r in 1:length(ratios)){
      W[,(d[1]+1):ncol(W[,,1]),r] <- W[,(d[1]+1):ncol(W[,,1]),r] * ratios[r]
    }
    
  }

 return(W)
    
}

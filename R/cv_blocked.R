#' @export
cv_blocked <- function(B, Z, Y, W, Ak, k, d, lambda1, t1, t2, eps,intercept=FALSE, cv, nfolds){
  
  # B  <- object@B
  # Z  <- t(as.matrix(object@A))
  # Y  <- t(as.matrix(object@b))
  # W  <- object@W
  # Ak  <- object@Ak
  # k  <- object@k
  # d  <- object@d
  # lambda1  <- object@lambda1
  # t1  <- object@t1
  # t2  <- object@t2
  # eps  <- 1e-3
  # intercept  <- object@intercept
  # cv  <- object@cv
  # nfolds  <- object@nfolds
  # B = B_com
  # Z = t(Z_com)
  # Y = t(object@b)
  # W = W_com
  
  make_folds  <- function(x,nfolds) split(x, cut(seq_along(x), nfolds, labels = FALSE))
  final_tmpt <- cumsum(unlist(lapply(Ak,function(x){nrow(x)})))
  first_tmpt <- c(1,(final_tmpt[-length(final_tmpt)]+1))
  subj_indx_list <- lapply(seq_along(first_tmpt), function(g){first_tmpt[g]:final_tmpt[g]})
  cv_list <- lapply(subj_indx_list,function(g){make_folds(g,nfolds)})
  #MSFE <- matrix(NA, nrow = nfolds, ncol = nrow(lambda1)*length(ratios))
  MSFE <- matrix(NA, nrow = nfolds, ncol = nrow(lambda1)*dim(W)[3])
  pb   <- txtProgressBar(1, nfolds, style=3)
  
  for(fold_id in 1:nfolds){ # fold_id <- 1
    
    setTxtProgressBar(pb, fold_id)
    
    test_idx  <- unlist(lapply(1:k,function(a){cv_list[[a]][fold_id]}))   
    train_idx <- unlist(lapply(1:k,function(a){cv_list[[a]][-fold_id]})) 
    
    # dim(B)
    # dim(Z[,train_idx])
    # dim(Y[,train_idx])
    # dim(W)
   
    # dataset1 dataset1          
    # 5       21      900 
    # [1]  20 267
    # [1]   5 267
    # [1]  5 20 30
    
    beta <- wlasso(B, Z[,train_idx], Y[,train_idx], W, k, d, lambda1,eps,intercept)
    
    # beta <- multivar:::wlasso(B[1,,,drop=F], Z[,train_idx], Y[1,train_idx,drop=F], W[1,,,drop=F], k, d, lambda1,eps,intercept)
    # beta <- multivar:::wlasso(B, Z[,train_idx], Y[1,train_idx,drop=F], W, k, d, lambda1,eps,intercept)
    
    # Calculate h-step MSFE for each penalty parameter
    for (ii in 1:dim(beta)[3]) {
      #MSFE[fold_id,ii] <- norm2(Y[,test_idx,drop=F]- beta[,-1,ii] %*% Z[,test_idx,drop=F] )^2
      MSFE[fold_id,ii] <- norm2(Y[,test_idx,drop=F]- beta[,,ii] %*% Z[,test_idx,drop=F] )^2
    }
  }
  
  beta <- wlasso(B, Z, Y, W, k, d, lambda1,eps,intercept)

  return(list(beta,MSFE))
  
  
}

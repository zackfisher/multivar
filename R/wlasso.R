#' @keywords internal
wlasso <- function(B, Z, Y, W = NULL, k, d, lambda1,eps,intercept=FALSE){
  
  # B <- B[1,,,drop=F]
  # Z <- Z[,train_idx]
  # Y <- Y[1,train_idx,drop=F]
  # W <- W[1,,,drop=F]
 
  if(is.null(W)){
    W <- array(1,dim=c(d[1],ncol(Z),1))
  }
  
  
  if(!"matrix"%in%class(Y)){
    Y <- matrix(Y,nrow=1)
  }

  Y     <- t(Y)

  tk <- 1/max(Mod(eigen(Z%*%t(Z),only.values=TRUE)$values))

  BFOO1 <- B[,,1] 

  beta <- lamloopFISTA(
    B, 
    Y, 
    Z, 
    W,
    as.matrix(lambda1),
    eps,
    #as.matrix(YMean), 
    #as.matrix(ZMean), 
    as.matrix(BFOO1),
    tk)

  return(beta)
}









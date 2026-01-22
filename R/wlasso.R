#' @export
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

  # begin intercept comment out
  #
  # if(intercept){
  #   
  #   YMean <- c(apply(Y, 1, mean))
  #   ZMean <- c(apply(Z, 1, mean))
  #   Y <- Y - YMean %*% t(c(rep(1, ncol(Y))))
  #   Z <- Z - ZMean %*% t(c(rep(1, ncol(Z))))
  #   
  # }else{
  #   YMean <- rep(0,nrow(Y))
  #   ZMean <- rep(0,nrow(Z))
  # }
  #
  # nc <- apply(B,3,ncol)[1
  # BFOO1 <- B[,2:nc,1] 
  # BFOO  <- B[,2:nc,,drop=F] 
  #
  # end
  
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









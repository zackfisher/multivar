#' Simulate a stationary Vector Autoregressive (VAR) time series.
#'
#' @param n An integer giving the number of timepoints.
#' @param phi A d x d transition matrix.
#' @param sigma A d x d innovation covariance matrix.
#' @param bubrn Burn-in period to be discarded. Default is 500.
#' @param tvp A d x d logical matrix indicating if a dynamic is time varying.
#' @keywords var simulate internal
#' @noRd
var_sim = function(n, phi, sigma, burn = 500,tvp=NULL,intercept=0){
  k    <- dim(phi)[1]
  p    <- dim(phi)[2]/k
  inno <- MASS::mvrnorm(n=n+burn, rep(0, k), sigma)
  init <- MASS::mvrnorm(n=p, rep(0, k), sigma)
  init <- matrix(init, nrow=p)
	j    <- 1
	id   <- seq(from= j+p-1, to = j, by=-1)
  Y    <-  matrix(0, (n+burn), k)
  
  phi <- replicate(n+burn, phi, simplify=FALSE)
  
  if(!is.null(tvp)){
    for(j in 1:nrow(tvp)){
      for(m in 1:ncol(tvp)){
        if(tvp[j,m]){
          
          tvp_series <- cumsum(sample(c(-1, 1), (n+burn), TRUE))
          tvp_series <- smooth(scales::rescale(tvp_series, to = c(-1, 1)))
          tvp_series <- cbind(tvp_series, index=c(1:(n+burn)))
          #tvp_series <- loess(tvp_series ~ index, data = as.data.frame(tvp_series), span=0.6)$fitted
          
          phi <- lapply(1:length(phi), function(v){
            phi[[v]][j,m] <- tvp_series[v]
            phi[[v]]
          })
          
        }
      }
    }
  } 
  
  
  for(r in 1:(n+burn)){
  	Y[r,] <-  intercept + phi[[r]] %*% as.vector(t(init[id,])) + inno[r,]
    init <- rbind(init[-1,], Y[r,])
  }
  
  if(is.null(tvp)){
    return(Y[-(1:burn),])
  } else {
    dat <- list(
      data = Y[-(1:burn),],
      phi  = phi[-(1:burn)]
    )
    return(dat)
  }
  
  
}

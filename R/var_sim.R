#' Simulate a stationary Vector Autoregressive (VAR) time series.
#'
#' @param n An integer giving the number of timepoints.
#' @param phi A d x d transition matrix.
#' @param sigma A d x d innovation covariance matrix.
#' @param bubrn Burn-in period to be discarded. Default is 500.
#' @keywords var simulate
#'
#' @examples
#'
#' theta    <- diag(c(.7,.8,.9,.6,.7,.9))
#' data     <- var_sim(100, theta, diag(.1,6))
#' datalag  <- embed(data, 2)
#' b        <- datalag[,1:6]
#' A        <- datalag[,7:12]
#' A_est    <- fista_sparse(A, b, 1, theta, niter = 10, backtrack = TRUE)$out.x
#' var_forecast(t(b), 2, A_est)
#'
#' @export
var_sim = function(n, phi, sigma, burn = 500){
  k    <- dim(phi)[1]
  p    <- dim(phi)[2]/k
  inno <- MASS::mvrnorm(n=n+burn, rep(0, k), sigma)
  init <- MASS::mvrnorm(n=p, rep(0, k), sigma)
  init <- matrix(init, nrow=p)
	j    <- 1
	id   <- seq(from= j+p-1, to = j, by=-1)
  Y    <-  matrix(0, (n+burn), k)
  for(r in 1:(n+burn)){
  	Y[r,] <-  phi %*% as.vector(t(init[id,])) + inno[r,]
     init <- rbind(init[-1,], Y[r,])
  }
  return(Y[-(1:burn),])
}

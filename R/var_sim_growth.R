#' Simulate a nonstationary Vector Autoregressive (VAR) time series.
#'
#' @param n An integer giving the number of timepoints.
#' @param phi A list of length n of d x d transition matrix.
#' @param sigma A d x d innovation covariance matrix.
#' @param burn Burn-in period to be discarded. Default is 500.
#' @param intercept Intercept term. Can be:
#'   - Scalar (same for all variables and time)
#'   - Vector of length d (different per variable, constant over time)
#'   - List of length n containing vectors (time-varying intercepts)
#' @keywords var simulate internal
#' @noRd
var_sim_growth = function(n, phi, sigma, burn = 500, intercept = 0){
  k    <- dim(phi[[1]])[1]
  p    <- dim(phi[[1]])[2]/k
  inno <- MASS::mvrnorm(n=n+burn, rep(0, k), sigma)
  init <- MASS::mvrnorm(n=p, rep(0, k), sigma)
  init <- matrix(init, nrow=p)
	j    <- 1
	id   <- seq(from= j+p-1, to = j, by=-1)
  Y    <-  matrix(0, (n+burn), k)

  phi_burn <- replicate(burn, phi[[1]], simplify=FALSE)
  phi      <- append(phi_burn, phi)

  # Handle intercept
  # If intercept is a list (time-varying), expand for burn-in
  if(is.list(intercept)){
    intercept_burn <- replicate(burn, intercept[[1]], simplify=FALSE)
    intercept_full <- append(intercept_burn, intercept)
  } else {
    # Scalar or vector - constant over time
    intercept_full <- intercept
  }

  for(r in 1:(n+burn)){
    # Get intercept for this time point
    if(is.list(intercept_full)){
      int_r <- intercept_full[[r]]
    } else {
      int_r <- intercept_full
    }

  	Y[r,] <-  int_r + phi[[r]] %*% as.vector(t(init[id,])) + inno[r,]
    init  <- rbind(init[-1,], Y[r,])
  }

  return(Y[-(1:burn),])

}

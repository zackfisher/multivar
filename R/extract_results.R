#' @export
extract_results <- function(x){
  #x <- r
  MSFE_mean    <- colMeans(x$MSFE)
  msfe_min_idx <- which.min(MSFE_mean)
  B <- x$beta[,,msfe_min_idx]

  mats <- breakup_transition(B, x$obj@Ak, x$obj@ndk, x$obj@intercept)
  
  return(mats)
}

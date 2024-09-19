#' Simulate a time-varying VAR model.
#'
#' @param n Integer. The time series length. 
#' @param sigma Matrix. Innovation covariance matrix.
#' @param mat_total Matrix. Transition matrix. 
#' @param mat_tvp A logical matrix indicating which parameters are time-varying.
#' @keywords var time-varying simulate
#' @examples
#' @export
tvp_sim <- function(
  n, 
  sigma, 
  mat_total,
  mat_tvp,
  int=0){
  
  dat <- lapply(mat_total, function(x) {var_sim(n, x, sigma,tvp=mat_tvp,int=int)})
  
  data <- lapply(lapply(dat,"[[", "data"), function(df){colnames(df) <- paste0("V",1:ncol(df)); df})
    
  mat_list <- list(
    mat_ind_final  = lapply(dat,"[[", "phi"),
    data = data
  )
    
  return(mat_list)
  
}


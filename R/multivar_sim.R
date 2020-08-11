#' Simulate multivar data.
#'
#' @param k Integer. The number of individuals (or datsets) to be generated.
#' @param d Integer. The number of variables per dataset. For now this will be constant across individuals. 
#' @param n Integer. The time series length. 
#' @param prop_fill_com Numeric. The proportion of nonzero paths in the common transition matrix.
#' @param prop_fill_ind Numeric. The proportion of nonzero unique (not in the common transition matrix or transition matrix of other individuals) paths in each individual transition matrix.
#' @param ub Numeric. The lower bound for individual elements of the transition matrices.
#' @param lb Numeric. The upper bound for individual elements of the transition matrices.
#' @param sigma Matrix. The (population) innovation covariance matrix.
#'
#' @keywords var multivar simulate
#'
#' @examples
#'
#' k <- 3
#' d <- 5
#' n <- 100
#' prop_fill_com <- .2
#' prop_fill_ind <- .2
#' lb <- 0.1
#' ub <- 0.7
#' sigma <- diag(0.1,d)
#' data <- multivar_sim(k, d, n, prop_fill_com, prop_fill_ind, lb, ub,sigma)$data
#'
#' @export
multivar_sim <- function(k, d, n, prop_fill_com, prop_fill_ind, lb, ub, sigma){
  
  densities             <- c(prop_fill_com, rep(prop_fill_ind,k))        
  n_elem_prop_fill_com  <- round(d^2*prop_fill_com) 
  n_elem_prop_fill_ind  <- round(d^2*prop_fill_ind) 
  n_elem_total_all_subj <- n_elem_prop_fill_ind*k
  n_elem_total_all      <- n_elem_total_all_subj + n_elem_prop_fill_com
  
  if(n_elem_total_all_subj > (d^2 - n_elem_prop_fill_com)){
    stop(paste0("multivar ERROR: The arguments prop_fill_com and prop_fill_ind are not well specified."))
  }
  
  permuted_elements <- sample(c(1:d^2))
  idx_grp           <- permuted_elements[1:n_elem_prop_fill_com]
  sub_elements      <- permuted_elements[(n_elem_prop_fill_com+1):(n_elem_prop_fill_com+n_elem_total_all_subj)]
  idx_sub           <- split(sub_elements, ceiling(seq_along(sub_elements)/n_elem_prop_fill_ind))
  idx_all           <- append(list(idx_grp),unname(idx_sub))
  max_eigs          <- TRUE
  
  while(max_eigs){
    
    mats <- replicate(k+1, matrix(0,d,d), simplify = FALSE)
    mats[[1]][idx_all[[1]]]<- d.fun(length(idx_all[[1]]))
      
    mats <- lapply(seq_along(mats), function(i){
      if(i == 1){
        mats[[i]]
      } else {
        mats[[i]][idx_all[[i]]] <- runif(length(idx_all[[i]]), lb, ub)
        mats[[i]] <- mats[[i]] + mats[[1]]
        mats[[i]]
      }
      
    })
    
    max_eigs <- any(unlist(lapply(mats[-1], function(x){max(abs(eigen(x)$values))})) > .9)
    
  }
  
  data <- lapply(mats[-1], function(x) {var_sim(n, x, sigma)})
  
  mat_list <- list(
    mat_com        = mats[[1]],
    mat_ind_unique = lapply(1:k, function(i){mats[[i+1]] - mats[[1]]}),
    mat_ind_final  = mats[-1],
    data = data
  )
  
  return(mat_list)
  
}


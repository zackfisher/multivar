#' Simulate multivar data.
#'
#' @param k Integer. The number of individuals (or datasets) to be generated.
#' @param d Integer. The number of variables per dataset. For now this will be constant across individuals. 
#' @param n Integer. The time series length. 
#' @param prop_fill_com Numeric. The proportion of nonzero paths in the common transition matrix.
#' @param prop_fill_ind Numeric. The proportion of nonzero unique (not in the common transition matrix or transition matrix of other individuals) paths in each individual transition matrix.
#' @param ub Numeric. The lower bound for individual elements of the transition matrices.
#' @param lb Numeric. The upper bound for individual elements of the transition matrices.
#' @param sigma Matrix. The (population) innovation covariance matrix.
#' @param unique_overlap Logical. Default is FALSE. Whether the unique portion should be completely unique (no overlap) or randomly chosen.
#' @param mat_common Matrix. A common effects transition matrix (if known).
#' @param mat_unique List. A list of unique effects transition matrix (if known).
#' @param mat_total List. A list of total effects transition matrix (if known).
#' @param diag Logical. Default is FALSE. Should diagonal elements be filled first for common elements.
#' @param intercept List. Default is NULL. A list of length K containing numeric vectors of length d representing the intercept values.
#' @keywords var multivar simulate
#' @examples
#' k <- 3
#' d <- 10
#' n <- 20
#' prop_fill_com <- .1
#' prop_fill_ind <- .05
#' lb <- 0.1
#' ub <- 0.5
#' sigma <- diag(d)
#' data <- multivar_sim(k, d, n, prop_fill_com, prop_fill_ind, lb, ub,sigma)$data
#' @export
multivar_sim <- function(
    k, 
    d, 
    n, 
    prop_fill_com, 
    prop_fill_ind, 
    lb, 
    ub,
    sigma, 
    unique_overlap = FALSE, 
    mat_common = NULL, 
    mat_unique = NULL, 
    mat_total = NULL, 
    diag = FALSE,
    intercept = NULL){
  
  if(is.null(intercept)){
    intercept <- replicate(k ,rep(0, d), simplify = FALSE)
  }
  
  # If the user did NOT supply a common matrix AND did NOT supply total matrices,
  # we will randomly generate the common + individual transition matrices here.
  if (is.null(mat_common) & is.null(mat_total)){
    
    # Helper function: draw n random coefficients uniformly from [lb, ub].
    d.fun                 <- function(n) runif(n,lb,ub)
    
    # Vector of sparsity “densities” (common first, then each individual).
    # (Note: not used later; kept as a bookkeeping object.)
    densities             <- c(prop_fill_com, rep(prop_fill_ind,k))        
    
    # Convert proportions into counts of nonzero elements (out of d^2 entries).
    n_elem_prop_fill_com  <- round(d^2*prop_fill_com) 
    n_elem_prop_fill_ind  <- round(d^2*prop_fill_ind) 
    
    # Total number of unique nonzero entries across all subjects (if no overlap).
    n_elem_total_all_subj <- n_elem_prop_fill_ind*k
    
    # Total number of nonzero “slots” needed = common + all subject-unique.
    n_elem_total_all      <- n_elem_total_all_subj + n_elem_prop_fill_com
    
    # If unique_overlap is FALSE, then the unique entries for all subjects must fit
    # into the complement of the common entries (i.e., no overlaps allowed).
    if(!unique_overlap & n_elem_total_all_subj > (d^2 - n_elem_prop_fill_com)){
      stop(paste0("multivar ERROR: The arguments prop_fill_com and prop_fill_ind are not well specified."))
    }
    
    # Flag controlling a rejection-sampling loop to enforce stability/stationarity.
    # We will keep drawing matrices until every subject matrix has max |eigenvalue| <= .95.
    max_eigs <- TRUE
    
    # Keep generating matrices until the spectral radius constraint is satisfied.
    while(max_eigs){
      # Optionally prioritize diagonal positions when choosing nonzero indices.
      if(diag){
        # Create a d-by-d matrix of indices 1..d^2, then vectorize it.
        vec_A             <- c(matrix(c(1:d^2),d,d))    # indices for elements of A
        
        # Extract the diagonal indices (in that same indexing scheme).
        vec_A_diag        <- diag(matrix(c(1:d^2),d,d)) # indices for diagonal elements of A
        
        # All remaining (non-diagonal) indices.
        vec_A_nondiag     <- setdiff(vec_A,vec_A_diag)  # indices for non-diagonal elements of A
        
        # Build a permutation that lists all diagonal indices first, then a random order of the rest.
        permuted_elements <- c(vec_A_diag,sample(vec_A_nondiag))
      } else {
        # Otherwise, just randomly permute all d^2 indices.
        permuted_elements <- sample(c(1:d^2))
      }
      
      # CASE 1: no common effects (prop_fill_com == 0), but there ARE individual unique effects.
      if(n_elem_prop_fill_com == 0 & n_elem_prop_fill_ind != 0){
        
        # Group/common index set is empty (represented as 0 here).
        idx_grp <- 0
        
        # If we require unique entries across subjects (no overlap), allocate disjoint blocks
        # of indices to subjects by slicing the permuted list.
        if(!unique_overlap){
          sub_elements <- permuted_elements[(1):(n_elem_total_all_subj)]
          idx_sub <- split(sub_elements, ceiling(seq_along(sub_elements)/n_elem_prop_fill_ind))
        } else {
          # If overlap is allowed, each subject samples its own unique indices (with possible overlap).
          sub_elements <- permuted_elements[(1):length(permuted_elements)]
          idx_sub <- lapply(1:k, function(v){sample(sub_elements,size = n_elem_prop_fill_ind)})
        }
        
        # Combine group indices + subject indices into one list (length k+1).
        idx_all <- append(list(idx_grp),unname(idx_sub))
        
        # Initialize k+1 matrices (first = “common”; the rest = subject totals) as all zeros.
        mats <- replicate(k+1, matrix(0,d,d), simplify = FALSE)
        
        # Fill each subject matrix at its chosen indices with random coefficients in [lb, ub].
        # The first matrix (common) stays as-is (all zeros) in this case.
        mats <- lapply(seq_along(mats), function(i){
          if(i == 1){
            mats[[i]]
          } else {
            mats[[i]][idx_all[[i]]] <- runif(length(idx_all[[i]]), lb, ub)
            mats[[i]]
          }
        })
        
        # CASE 2: there ARE common effects, but no individual unique effects (prop_fill_ind == 0).
      } else if (n_elem_prop_fill_com != 0 & n_elem_prop_fill_ind == 0){
        
        # Select the first n_elem_prop_fill_com indices for the common matrix.
        idx_grp <- permuted_elements[1:n_elem_prop_fill_com]
        
        # Store only the common index set in a list.
        idx_all <- list(idx_grp)
        
        # Initialize k+1 matrices as zeros.
        mats <- replicate(k+1, matrix(0,d,d), simplify = FALSE)
        
        # Fill the common matrix’s chosen indices with random coefficients.
        mats[[1]][idx_all[[1]]]<- d.fun(length(idx_all[[1]]))
        
        # Copy the common matrix into each subject’s total matrix (so all subjects identical).
        mats <- lapply(seq_along(mats), function(i){
          if(i == 1){
            mats[[i]]
          } else {
            mats[[i]] <-  mats[[1]]
            mats[[i]]
          }
        })
        
        # CASE 3: there ARE common effects AND there ARE individual unique effects.
      } else if (n_elem_prop_fill_com != 0 & n_elem_prop_fill_ind != 0){
        
        # Select indices for the common matrix.
        idx_grp <- permuted_elements[1:n_elem_prop_fill_com]
        
        # Select indices for the individual-unique parts.
        if(!unique_overlap){
          # If no overlap allowed, take a disjoint slice after the common indices.
          sub_elements <- permuted_elements[(n_elem_prop_fill_com+1):(n_elem_prop_fill_com+n_elem_total_all_subj)]
          idx_sub <- split(sub_elements, ceiling(seq_along(sub_elements)/n_elem_prop_fill_ind))
        } else {
          # If overlap allowed, each subject samples from the remaining indices (after common).
          sub_elements <- permuted_elements[(n_elem_prop_fill_com+1):length(permuted_elements)]
          idx_sub <- lapply(1:k, function(v){sample(sub_elements,size = n_elem_prop_fill_ind)})
        }
        
        # Combine common indices + subject-unique indices into one list (length k+1).
        idx_all <- append(list(idx_grp),unname(idx_sub))
        
        # Initialize matrices (common + k subject totals) as zeros.
        mats <- replicate(k+1, matrix(0,d,d), simplify = FALSE)
        
        # Fill the common matrix at its indices.
        mats[[1]][idx_all[[1]]]<- d.fun(length(idx_all[[1]]))
        
        # For each subject: fill unique indices, then add the common matrix to get the total.
        mats <- lapply(seq_along(mats), function(i){
          if(i == 1){
            mats[[i]]
          } else {
            mats[[i]][idx_all[[i]]] <- runif(length(idx_all[[i]]), lb, ub)
            mats[[i]] <- mats[[i]] + mats[[1]]
            mats[[i]]
          }
        })
        
      } else {
        # If both proportions are zero, there is nothing to simulate.
        stop(paste0("multivar ERROR: One of the arguments prop_fill_com and prop_fill_ind must be nonzero."))
      }
      
      # Stability check: compute max absolute eigenvalue (“spectral radius”) for each subject matrix,
      # and repeat if any subject’s max |eigenvalue| exceeds 0.95.
      max_eigs <- any(unlist(lapply(mats[-1], function(x){max(abs(eigen(x)$values))})) > .95)
    }
    
    # Simulate VAR time series data for each subject using its transition matrix and innovation covariance.
    data <- lapply(seq_along(mats[-1]), function(x_i) {var_sim(n, mats[-1][[x_i]], sigma, intercept=intercept[[x_i]])})
    
    # Assign default column names V1, V2, ..., Vd for each subject’s dataset.
    data <- lapply(data, function(df){colnames(df) <- paste0("V",1:ncol(df)); df})
    
    # Package outputs: common matrix, unique (subject - common), total (subject) matrices, and data.
    mat_list <- list(
      mat_com        = mats[[1]],
      mat_ind_unique = lapply(1:k, function(i){mats[[i+1]] - mats[[1]]}),
      mat_ind_final  = mats[-1],
      data = data
    )
    
    
  } else {
    
    # If the user supplied some matrices, we use them instead of generating new ones.
    # If mat_total is missing but mat_common + mat_unique are provided, build mat_total.
    if(is.null(mat_total)){
      mat_total <- lapply(1:length(mat_unique), function(i){mat_common + mat_unique[[i]]})
    }
    
    # Simulate VAR time series data for each provided total transition matrix.
    data <- lapply(seq_along(mat_total), function(x_i) {var_sim(n, mat_total[[x_i]], sigma, intercept = intercept[[x_i]])})
    
    # Assign default column names V1, V2, ..., Vd for each dataset.
    data <- lapply(data, function(df){colnames(df) <- paste0("V",1:ncol(df)); df})
    
    # Package outputs consistently with the generated-matrix branch.
    mat_list <- list(
      mat_com        = mat_common,
      mat_ind_unique = mat_unique,
      mat_ind_final  = mat_total,
      data = data
    )
    
    
  }
  
  # Return the list containing matrices + simulated datasets.
  return(mat_list)
  
}





#' Simulate multivar data.
#'
#' @param k Integer. The number of individuals (or datasets) to be generated.
#' @param d Integer. The number of variables per dataset. For now this will be constant across individuals. 
#' @param n Integer. The time series length. 
#' @param prop_fill_com Numeric. The proportion of nonzero paths in the common transition matrix.
#' @param prop_fill_ind Numeric. The proportion of nonzero unique (not in the common transition matrix or transition matrix of other individuals) paths in each individual transition matrix.
#' @param ub Numeric. The lower bound for individual elements of the transition matrices.
#' @param lb Numeric. The upper bound for individual elements of the transition matrices.
#' @param sigma Matrix. The (population) innovation covariance matrix.
#' @param unique_overlap Logical. Default is FALSE. Whether the unique portion should be completely unique (no overlap) or randomly chosen.
#' @param mat_common Matrix. A common effects transition matrix (if known).
#' @param mat_unique List. A list of unique effects transition matrix (if known).
#' @param mat_total List. A list of total effects transition matrix (if known).
#' @param diag Logical. Default is FALSE. Should diagonal elements be filled first for common elements.
#' @keywords var multivar simulate
#' @examples
#' k <- 3
#' d <- 10
#' n <- 20
#' prop_fill_com <- .1
#' prop_fill_ind <- .05
#' lb <- 0.1
#' ub <- 0.5
#' sigma <- diag(d)
#' data <- multivar_sim(k, d, n, prop_fill_com, prop_fill_ind, lb, ub,sigma)$data
#' @export
multivar_sim_old <- function(
  k, 
  d, 
  n, 
  prop_fill_com, 
  prop_fill_ind, 
  lb, 
  ub,
  sigma, 
  unique_overlap = FALSE, 
  mat_common = NULL, 
  mat_unique = NULL, 
  mat_total = NULL, 
  diag = FALSE){
  

  if (is.null(mat_common) & is.null(mat_total)){
  
    d.fun                 <- function(n) runif(n,lb,ub)
    densities             <- c(prop_fill_com, rep(prop_fill_ind,k))        
    n_elem_prop_fill_com  <- round(d^2*prop_fill_com) 
    n_elem_prop_fill_ind  <- round(d^2*prop_fill_ind) 
    n_elem_total_all_subj <- n_elem_prop_fill_ind*k
    n_elem_total_all      <- n_elem_total_all_subj + n_elem_prop_fill_com
    
    if(!unique_overlap & n_elem_total_all_subj > (d^2 - n_elem_prop_fill_com)){
      stop(paste0("multivar ERROR: The arguments prop_fill_com and prop_fill_ind are not well specified."))
    }
    
    max_eigs <- TRUE
    
    while(max_eigs){
      if(diag){
        vec_A             <- c(matrix(c(1:d^2),d,d))    # indices for elements of A
        vec_A_diag        <- diag(matrix(c(1:d^2),d,d)) # indices for diagonal elements of A
        vec_A_nondiag     <- setdiff(vec_A,vec_A_diag)  # indices for non-diagonal elements of A
        permuted_elements <- c(vec_A_diag,sample(vec_A_nondiag))
      } else {
        permuted_elements <- sample(c(1:d^2))
      }
      
      if(n_elem_prop_fill_com == 0 & n_elem_prop_fill_ind != 0){
        
        idx_grp <- 0
        if(!unique_overlap){
          sub_elements <- permuted_elements[(1):(n_elem_total_all_subj)]
          idx_sub <- split(sub_elements, ceiling(seq_along(sub_elements)/n_elem_prop_fill_ind))
        } else {
          sub_elements <- permuted_elements[(1):length(permuted_elements)]
          idx_sub <- lapply(1:k, function(v){sample(sub_elements,size = n_elem_prop_fill_ind)})
        }
        idx_all <- append(list(idx_grp),unname(idx_sub))
        mats <- replicate(k+1, matrix(0,d,d), simplify = FALSE)
        mats <- lapply(seq_along(mats), function(i){
          if(i == 1){
            mats[[i]]
          } else {
            mats[[i]][idx_all[[i]]] <- runif(length(idx_all[[i]]), lb, ub)
            mats[[i]]
          }
        })
        
      } else if (n_elem_prop_fill_com != 0 & n_elem_prop_fill_ind == 0){
        
        idx_grp <- permuted_elements[1:n_elem_prop_fill_com]
        idx_all <- list(idx_grp)
        mats <- replicate(k+1, matrix(0,d,d), simplify = FALSE)
        mats[[1]][idx_all[[1]]]<- d.fun(length(idx_all[[1]]))
        mats <- lapply(seq_along(mats), function(i){
          if(i == 1){
            mats[[i]]
          } else {
            mats[[i]] <-  mats[[1]]
            mats[[i]]
          }
        })
        
      } else if (n_elem_prop_fill_com != 0 & n_elem_prop_fill_ind != 0){
        
        idx_grp <- permuted_elements[1:n_elem_prop_fill_com]
        if(!unique_overlap){
          sub_elements <- permuted_elements[(n_elem_prop_fill_com+1):(n_elem_prop_fill_com+n_elem_total_all_subj)]
          idx_sub <- split(sub_elements, ceiling(seq_along(sub_elements)/n_elem_prop_fill_ind))
        } else {
          sub_elements <- permuted_elements[(n_elem_prop_fill_com+1):length(permuted_elements)]
          idx_sub <- lapply(1:k, function(v){sample(sub_elements,size = n_elem_prop_fill_ind)})
        }
        idx_all <- append(list(idx_grp),unname(idx_sub))
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
        
      } else {
         stop(paste0("multivar ERROR: One of the arguments prop_fill_com and prop_fill_ind must be nonzero."))
      }

      max_eigs <- any(unlist(lapply(mats[-1], function(x){max(abs(eigen(x)$values))})) > .95)
    }
    
    data <- lapply(mats[-1], function(x) {var_sim(n, x, sigma)})
    data <- lapply(data, function(df){colnames(df) <- paste0("V",1:ncol(df)); df})
    
    mat_list <- list(
      mat_com        = mats[[1]],
      mat_ind_unique = lapply(1:k, function(i){mats[[i+1]] - mats[[1]]}),
      mat_ind_final  = mats[-1],
      data = data
    )
  
    
  } else {
    
    if(is.null(mat_total)){
      mat_total <- lapply(1:length(mat_unique), function(i){mat_common + mat_unique[[i]]})
    }
    
    data <- lapply(mat_total, function(x) {var_sim(n, x, sigma)})
    data <- lapply(data, function(df){colnames(df) <- paste0("V",1:ncol(df)); df})
    
    mat_list <- list(
      mat_com        = mat_common,
      mat_ind_unique = mat_unique,
      mat_ind_final  = mat_total,
      data = data
    )
    
    
  }
  
  return(mat_list)
  
}


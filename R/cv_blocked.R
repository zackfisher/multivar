#' Build Row Specification from Ak and Breaks
#'
#' Creates a minimal row specification when full spec is not available.
#' Used by legacy functions that don't have access to the full matrix_spec.
#'
#' @param Ak List of design matrices per subject
#' @param breaks List of period indices per subject (optional, for TVP models)
#'
#' @return List with rows$subject and optionally rows$period
#' @keywords internal
build_row_spec <- function(Ak, breaks = NULL) {
  k <- length(Ak)
  final_tmpt <- cumsum(unlist(lapply(Ak, function(x) nrow(x))))
  first_tmpt <- c(1, (final_tmpt[-length(final_tmpt)] + 1))

  # Build subject boundaries
  subject_rows <- lapply(seq_len(k), function(i) {
    c(first_tmpt[i], final_tmpt[i])
  })

  # Build period boundaries if breaks provided
  period_rows <- NULL
  if (!is.null(breaks)) {
    period_rows <- lapply(seq_len(k), function(g) {
      lapply(breaks[[g]], function(period_indices) {
        global_start <- period_indices[1] + first_tmpt[g] - 1
        global_end <- period_indices[length(period_indices)] + first_tmpt[g] - 1
        c(global_start, global_end)
      })
    })
  }

  list(rows = list(subject = subject_rows, period = period_rows))
}


#' Blocked Cross-Validation for multivar
#'
#' Performs k-fold blocked cross-validation for time series data.
#'
#' @param B Initial coefficient array
#' @param Z Design matrix (transposed, d x n)
#' @param Y Response matrix (transposed, d x n)
#' @param W Weight array for adaptive LASSO
#' @param Ak List of design matrices per subject (used to build row spec if spec not provided)
#' @param k Number of subjects
#' @param d Number of variables
#' @param lambda1 Lambda grid matrix
#' @param t1 Start indices (unused, kept for API compatibility)
#' @param t2 End indices (unused, kept for API compatibility)
#' @param eps Convergence tolerance
#' @param intercept Whether model includes intercepts
#' @param cv CV method (unused here, always "blocked")
#' @param nfolds Number of CV folds
#' @param tvp Whether time-varying parameters are used
#' @param breaks List of period indices per subject (required if tvp=TRUE)
#' @param spec Optional matrix_spec object. If not provided, a minimal spec is
#'             built from Ak and breaks.
#'
#' @return List with beta coefficients and MSFE matrix
#' @export
cv_blocked <- function(B, Z, Y, W, Ak, k, d, lambda1, t1, t2, eps,
                       intercept = FALSE, cv, nfolds, tvp = FALSE,
                       breaks = NULL, spec = NULL) {

  # Build spec if not provided (for legacy callers)
  if (is.null(spec)) {
    spec <- build_row_spec(Ak, breaks)
  }

  make_folds <- function(x, nfolds) split(x, cut(seq_along(x), nfolds, labels = FALSE))

  # Get subject row boundaries from spec

  first_tmpt <- sapply(seq_len(k), function(i) spec$rows$subject[[i]][1])
  final_tmpt <- sapply(seq_len(k), function(i) spec$rows$subject[[i]][2])

  # If TVP, create folds within each period separately to respect period boundaries
  if (tvp && !is.null(breaks)) {
    cv_list <- lapply(seq_len(k), function(g) {
      n_periods <- length(spec$rows$period[[g]])

      # Create folds within each period
      period_folds <- lapply(seq_len(n_periods), function(p) {
        period_range <- spec$rows$period[[g]][[p]]
        global_indices <- seq(period_range[1], period_range[2])
        make_folds(global_indices, nfolds)
      })

      # Combine folds across periods
      # For each fold_id, concatenate indices from the same fold across all periods
      combined_folds <- lapply(seq_len(nfolds), function(fold_id) {
        # Get fold_id from each period (may not exist if period is shorter than nfolds)
        fold_indices <- lapply(period_folds, function(pf) {
          if(fold_id <= length(pf)) {
            pf[[fold_id]]
          } else {
            integer(0)  # Empty if this fold doesn't exist for this period
          }
        })
        unlist(fold_indices)
      })

      combined_folds
    })
  } else {
    # Original behavior: create folds across entire time series
    subj_indx_list <- lapply(seq_along(first_tmpt), function(g){first_tmpt[g]:final_tmpt[g]})
    cv_list <- lapply(subj_indx_list,function(g){make_folds(g,nfolds)})
  }
  #MSFE <- matrix(NA, nrow = nfolds, ncol = nrow(lambda1)*length(ratios_unique))
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

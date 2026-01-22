## ------------------------------------------------------------
## cross-validation for N=1 adaptive VAR
## ------------------------------------------------------------
#' @export
cv.var <- function(object, lambda_best = "min", adaptive = TRUE, intercept = FALSE) {
  
  d      <- object@d[1L]
  nscen  <- length(object@ratios)
  
  ## ---- 3) (re)initialize coefs on residualized data ----
  if(adaptive){
    
    init_ind_dynamics <- estimate_initial_coefs(
      Ak        = object@Ak,
      bk        = object@bk,
      d         = object@d,
      k         = 1L,
      lassotype = object@lassotype,
      weightest = object@weightest,
      subgroup_membership = object@subgroup_membership,
      subgroup  = object@subgroup,
      nlambda1  = object@nlambda1,
      nlambda2  = object@nlambda2,
      tvp       = object@tvp,
      breaks    = object@breaks
    )$total_effects
    
  }
  
  object@ratios <- rep(1, length(object@ratios))
  
  ## Initial B slices per scenario (match Stage 1 convention)
  B_empty <- array(0, dim = c(d, d + 1L, dim(object@B)[3]))
  
  ## ---- 4) fit per group with ptLasso-style penalty factors encoded in Wk ----
  fit_k <- lapply(seq_len(object@k), function(k) {
    
    Wk <- array(0, dim = c(d, d, nscen))
    
    if(adaptive){
      
      adapower <- 1
      Ahat  <- as.matrix(init_ind_dynamics[[k]])   # d x (d+1)
      w_mat <- 1 / pmax(abs(Ahat)^adapower, .Machine$double.eps)  # d x d
      for (s in seq_len(nscen)) Wk[ , , s] <- w_mat
      
    } else {
      
      for (s in seq_len(nscen)) Wk[ , , s] <- 1
      
    }
    
    ## 4b) subject-specific lambda path on residualized (Y, X)
    lambda1k <- lambda_grid(
      depth     = object@depth,
      nlam      = object@nlambda1,
      Y         = t(as.matrix(object@bk[[k]])),     # d x T_k
      Z         = t(as.matrix(object@Ak[[k]])),# d x T_k
      W         = Wk,                          # d x d x nscen
      tol       = object@tol,
      intercept = object@intercept,
      lamadapt  = object@lamadapt,
      k         = 1L
    )
    
    ## 4c) blocked CV on the reduced per-group problem
    cvfit <- cv_multivar(
      B = B_empty,
      Z = t(as.matrix(object@Ak[[k]])),
      Y = t(as.matrix(object@bk[[k]])),
      W = Wk,
      Ak = list(as.matrix(object@Ak[[k]])),
      bk = list(as.matrix(object@bk[[k]])),
      1L,
      object@d,
      lambda1k,
      object@t1,
      object@t2,
      eps = object@tol,
      object@intercept,
      object@cv,
      object@nfolds
    )
    
    hyp      <- extract_multivar_hyperparams(object, cvfit)
    beta_arr <- cvfit[[1]]  # (d x (d+1) x nlam_all)
    MSFE     <- cvfit[[2]]  # (nfolds x nlam_all), flattened over scenarios
    
    #hyp[which(colMeans(MSFE)==min(colMeans(MSFE))),]
    
    lambda1_min <- hyp[which(colMeans(MSFE)==min(colMeans(MSFE))),"lambda1_value"][1]
    lambda1_min_msfe <- min(colMeans(MSFE))
    
    MSFE_min <- colMeans(MSFE)[which.min(colMeans(MSFE))]
    
    if(lambda_best == "min"){
      
      best_idx <- which.min(colMeans(MSFE))
      
    } else if (lambda_best == "1se"){
      
      cv_mean  <- colMeans(MSFE, na.rm=TRUE)
      cv_se    <- apply(MSFE, 2, sd, na.rm=TRUE) / sqrt(nrow(MSFE))
      imin     <- which.min(cv_mean)
      new_msfe <- cv_mean[imin] + cv_se[imin]
      lambda1_1se_msfe <- min(colMeans(MSFE)[colMeans(MSFE) > new_msfe])
      lambda1_1se <- hyp[hyp$MSFE == lambda1_1se_msfe,"lambda1_value"][1]
      best_idx <- which(colMeans(MSFE) == lambda1_1se_msfe)[1]
    }
    
    U_hat     <- beta_arr[ , , best_idx, drop = TRUE]         # d x (d+1)                        # d x (d+1)
    A_hat_k   <- U_hat                             # final subject-k matrix
    
    if(!intercept){
      A_hat_k <- A_hat_k[, -1, drop = FALSE]
    } 
    
    list(
      U           = U_hat,        # subject-specific refinement (intercept + slopes)
      A_hat       = A_hat_k,      # final A for subject k
      hyp         = hyp,
      MSFE        = MSFE,
      MSFE_min    = MSFE_min,
      best_idx    = best_idx,
      hyp         = hyp,
      method      = "pretrained_var_stage2"
    )
  })
  
  mats <- list(
    "common" = NULL,
    "subgrp" = NULL,
    "unique" = NULL,
    "total" = lapply(fit_k, function(x){ x$A_hat })
  )
  
  return(
    list(
      fit = fit_k,
      mats = mats
    )
  )
}

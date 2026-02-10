## ------------------------------------------------------------
## Stage 2: subject-specific fine-tuning with offset
## ------------------------------------------------------------
#' @keywords internal
pretrained_var_stage2 <- function(object, stage1, alpha = 0, stage2adaptive = FALSE, stage2ratios = FALSE, intercept = FALSE, lambda_best = "1se") {
  
  d      <- object@d[1L]
  nscen  <- length(object@ratios_unique)
  
  ## ---- 1) split Stage-1 into intercept and slopes (always use d x d for slopes) ----
  has_int_col <- (ncol(stage1$A_pre) == d + 1L)
  A_slopes    <- if (has_int_col) stage1$A_pre[, -1, drop = FALSE] else stage1$A_pre
  b0          <- if (has_int_col && isTRUE(object@intercept)) {
    stage1$A_pre[, 1, drop = FALSE]    # d x 1
  } else {
    matrix(0, nrow = d, ncol = 1)      # d x 1
  }
  
  ## ---- 2) residualize outcomes per group: Y_k - (1 - alpha) * (b0 + X_k A_slopes^T) ----
  bk_r <- lapply(seq_len(object@k), function(k) {
    Ak   <- object@Ak[[k]]                      # T_k x d
    Yk   <- object@bk[[k]]                      # T_k x d
    T_k  <- nrow(Ak)
    
    # replicate intercepts b0 across time: T_k x d (one per response)
    b1   <- matrix(rep(drop(b0), each = T_k), nrow = T_k, ncol = d, byrow = FALSE)
    
    off_k <- (1 - alpha) * (b1 + Ak %*% t(A_slopes))  # T_k x d
    Yk - off_k
  })
  
  ## ---- 3) (re)initialize coefs on residualized data ----
  if(stage2adaptive){
    
    init_ind_dynamics <- estimate_initial_coefs(
      Ak        = object@Ak,
      bk        = bk_r,
      d         = object@d,
      k         = 1L,
      lassotype = object@lassotype,
      weightest = object@weightest,
      subgroup_membership = object@subgroup_membership,
      subgroup  = object@subgroup,
      nlambda1  = object@nlambda1,
      tvp       = object@tvp,
      breaks    = object@breaks,
      intercept = object@intercept,
      nfolds    = object@nfolds
    )$total_effects
    
  }
  
  if(!stage2ratios){
    object@ratios_unique <- rep(1, length(object@ratios_unique))
  }
  
  ## Initial B slices per scenario (match Stage 1 convention)
  B_empty <- array(0, dim = c(d, d + 1L, dim(object@B)[3]))
  
  ## ---- 4) fit per group with ptLasso-style penalty factors encoded in Wk ----
  fit_k <- lapply(seq_len(object@k), function(k) {
    
    ## 4a) build penalty-factor matrix on slopes: 1 on Stage-1 support, 1/alpha off-support
    supp      <- (abs(A_slopes) > 0)                   # d x d
    alpha_eff <- if (alpha > 0) alpha else 1e-9        # stabilize alpha = 0
    pf_mat    <- matrix(1 / alpha_eff, nrow = d, ncol = d)
    pf_mat[supp] <- 1
    
    Wk <- array(0, dim = c(d, d, nscen))
    
    if(stage2adaptive){
      
      adapower <- 1
      Ahat  <- as.matrix(init_ind_dynamics[[k]])   # d x (d+1)
      w_mat <- 1 / pmax(abs(Ahat)^adapower, .Machine$double.eps)  # d x d
      for (s in seq_len(nscen)) Wk[ , , s] <- pf_mat * w_mat
      
    } else {
      
      for (s in seq_len(nscen)) Wk[ , , s] <- pf_mat  # * object@ratios_unique[s]
      
    }
    
    ## 4b) subject-specific lambda path on residualized (Y, X)
    lambda1k <- lambda_grid(
      depth     = object@depth,
      nlam      = object@nlambda1,
      Y         = t(as.matrix(bk_r[[k]])),     # d x T_k
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
      Y = t(as.matrix(bk_r[[k]])),
      W = Wk,
      Ak = list(as.matrix(object@Ak[[k]])),
      bk = list(as.matrix(bk_r[[k]])),
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
    
    U_hat     <- beta_arr[ , , best_idx, drop = TRUE]         # d x (d+1)
    A_common  <- cbind(b0, A_slopes)                          # d x (d+1)
    A_hat_k   <- A_common + U_hat                             # final subject-k matrix
    
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
  
  return(fit_k)
}

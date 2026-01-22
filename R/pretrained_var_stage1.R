#' Stage 1 for pretrained VAR: common block only
#'
#' Runs blocked CV on the reduced problem:
#'   Y_t = X_t^{(com)} A_0^T + e_t
#' where X_t^{(com)} are the first d columns of object@A.
#'
#' @export
pretrained_var_stage1 <- function(object, stage1adaptive = FALSE, stage1ratios = FALSE, lambda_best = "1se") {
  
  # notes on dimension:
  #   - W has length(ratios) slices
  #   - B has length(ratios)*length(lambda1) slices
  
  d      <- object@d[1L]
  #object@ratios <- rep(1, length(object@ratios))  # ignore ratios in stage 1
  nscen  <- length(object@ratios)
  
  # Common predictors (N x d)
  A_com <- object@A[, seq_len(d), drop = FALSE]
  
  # B for the reduced problem (d outcomes x [intercept + d slopes] x nscen)
  B_com <- object@B[, c(1:(d+1)), ,drop = FALSE]
  
  # ---------- weights ----------
  # Standard: uniform weights over slopes (intercept unpenalized by design)
  # Adaptive (optional): build d x d matrix and replicate across scenarios
  if (stage1adaptive) {
    
    adapower <- 1
    
    initcoefs <- estimate_initial_coefs(
      Ak        = list(do.call(rbind, object@Ak)),
      bk        = list(do.call(rbind, object@bk)),
      d         = object@d,
      k         = 1L,
      lassotype = object@lassotype,   # init with standard lasso
      weightest = object@weightest,
      subgroup_membership = object@subgroup_membership,
      subgroup  = object@subgroup,
      nlambda1  = object@nlambda1,
      nlambda2  = object@nlambda2,
      tvp       = object@tvp,
      breaks    = object@breaks
    )
    
    Ahat  <- as.matrix(initcoefs$total_effects[[1]])   # d x (d+1)
    w_mat <- 1 / pmax(abs(Ahat)^adapower, .Machine$double.eps)  # d x d
    
  } else {
    
    w_mat <- matrix(1, nrow = d, ncol = d)  # standard weights
    
  }
  
  W <- array(0, dim = c(d, d, nscen))
  
  if(stage1ratios){
    
    for (s in seq_len(nscen)) W[, , s] <- w_mat * object@ratios[s]
    
  } else {
    
    object@ratios <- rep(1, length(object@ratios))
    for (s in seq_len(nscen)) W[, , s] <- w_mat * 1
    
  }
 
  
  # Respect pendiag = FALSE by (nearly) unpenalizing diagonal
  #if (isFALSE(object@pendiag)) diag(w_mat) <- 1e-10

  # ---------- lambda grid on reduced problem ----------
  object@lambda1 <- lambda_grid(
    depth     = object@depth,
    nlam      = object@nlambda1,
    Y         = t(object@b),
    Z         = t(A_com),
    W         = W,
    tol       = object@tol,
    intercept = object@intercept,
    lamadapt  = object@lamadapt,
    k         = 1L
  )
  
  # ---------- blocked CV on reduced problem ----------
  fit <- cv_multivar(
    B = B_com,
    Z = as.matrix(t(A_com)),
    Y = as.matrix(t(object@b)),
    W = W,
    object@Ak,
    object@bk,
    k = 1L,
    object@d,
    object@lambda1,
    object@t1,
    object@t2,
    eps = object@tol,
    object@intercept,
    object@cv,
    object@nfolds
  )
  
  
  
  hyp      <- extract_multivar_hyperparams(object, fit)
  beta_arr <- fit[[1]]  # (d x (d+1) x nlam_all)
  MSFE     <- fit[[2]]  # (nfolds x nlam_all), flattened over scenarios
  
  #hyp[which(colMeans(MSFE)==min(colMeans(MSFE))),]
      
  lambda1_min <- hyp[which(colMeans(MSFE)==min(colMeans(MSFE))),"lambda1_value"][1]
  lambda1_min_msfe <- min(colMeans(MSFE))
  MSFE_min <- colMeans(MSFE)[which.min(colMeans(MSFE))]
  
  cv_mean  <- colMeans(MSFE, na.rm=TRUE)
  cv_se    <- apply(MSFE, 2, sd, na.rm=TRUE) / sqrt(nrow(MSFE))
  imin     <- which.min(cv_mean)
  new_msfe <- cv_mean[imin] + cv_se[imin]
  
  # find the closest element of colMeans(MSFE) to new_msfe, without going below
  lambda1_1se_msfe <- min(colMeans(MSFE)[colMeans(MSFE) > new_msfe])
  lambda1_1se <- hyp[hyp$MSFE == lambda1_1se_msfe,"lambda1_value"][1]
  
 
  if(lambda_best == "min"){
    best_idx <- which.min(colMeans(MSFE))
  } else if (lambda_best == "1se"){

    # there are multiple in this case but we just need one
    # to see this:  
    # beta_arr[, , (colMeans(MSFE) == lambda1_1se_msfe), drop = TRUE] 
    best_idx <- which(colMeans(MSFE) == lambda1_1se_msfe)[1]
  }
  
  A_pre    <- beta_arr[, , best_idx, drop = TRUE]   # d x (d+1)
  
  # for(jj in 1:dim(beta_arr)[3]){
  #   print(paste0(jj,":",sum(beta_arr[,,jj]!=0)))
  # }
  # 
  # for(jj in 1:dim(beta_arr)[3]){
  #   print(paste0(jj,":",mean(MSFE[,jj])))
  # }
  
  #find all the duplciated matrices in beta_arr[,,]
  # unique_matrices <- list()
  # for(jj in 1:dim(beta_arr)[3]){
  #   mat_jj <- beta_arr[,,jj]
  #   is_duplicate <- FALSE
  #   for(um in unique_matrices){
  #     if(all(mat_jj == um)){
  #       is_duplicate <- TRUE
  #       break
  #     }
  #   }
  #   if(!is_duplicate){
  #     unique_matrices[[length(unique_matrices)+1]] <- mat_jj
  #   }
  # }
  #print(paste0("Number of unique matrices in stage 1: ", length(unique_matrices)))
  
  # check number of unique rows in MSFE
  #unique_msfe <- unique(colMeans(MSFE))
  #print(paste0("Number of unique MSFE values in stage 1: ", length(unique_msfe)))
  
  out <- list(
    A_pre    = A_pre,
    A_com    = A_com,
    hyp      = hyp,
    MSFE     = MSFE,
    MSFE_min = MSFE_min,
    best_idx = best_idx,
    obj      = object,
    method   = "pretrained_var_stage1"
  )
  class(out) <- c("multivar.pretrain.stage1", class(out))
  out
}




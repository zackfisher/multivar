#' Stage 1 (glmnet): common block only, with optional per-equation intercepts
#'
#' Fits the common block using cv.glmnet on a block-diagonal design.
#' Returns A_pre (d x (d+1)); if include_intercept = FALSE, the first column is 0.
#'
#' @param object multivar object with slots: A (N x p), b (N x d), d (length-1 or length-d),
#'               nfolds, etc.
#' @param include_intercept logical; if TRUE, add one unpenalized intercept per outcome (default TRUE).
#' @param standardize logical; passed to glmnet::cv.glmnet (default FALSE to mirror multivar defaults).
#' @param foldid optional CV fold id. If length N, it will be replicated for each outcome (length N*d).
#'               If length N*d, used as-is. Otherwise ignored (glmnet will make its own).
#' @param lambda optional numeric vector of lambdas to force glmnet to use (helps match multivar paths).
#' @return list(A_pre, lambda_best, cvfit, include_intercept)
#' @importFrom glmnet cv.glmnet
#' @export
pretrained_var_stage1_glmnet <- function(object,
                                         include_intercept = FALSE,
                                         standardize = FALSE,
                                         foldid = NULL,
                                         lambda = NULL,
                                         lambest = "min") {
  
  # dims
  d  <- object@d[1L]
  A_com <- object@A[, seq_len(d), drop = FALSE]  # N x d
  Y     <- object@b                              # N x d
  N     <- nrow(A_com)
  
  # stack y so rows 1..N are outcome 1, then outcome 2, ..., outcome d
  y_stack <- as.vector(as.matrix(Y))  # length N*d (column-major)
  
  # block-diagonal design: kronecker(I_d, A_com) -> (N*d) x (d*d)
  X_blk <- kronecker(diag(d), A_com)
  
  if (include_intercept) {
    # per-outcome intercept columns: (N*d) x d
    Inter_blk <- kronecker(diag(d), matrix(1, nrow = N, ncol = 1))
    X <- cbind(Inter_blk, X_blk)
    
    # unpenalized per-outcome intercepts, unit-penalized slopes
    pf <- c(rep(0, d), rep(1, d * d))
  } else {
    X  <- X_blk
    pf <- rep(1, d * d)
  }
  
  # fold id handling (optional): allow N or N*d length
  foldid_glm <- NULL
  if (!is.null(foldid)) {
    if (length(foldid) == N) {
      foldid_glm <- rep(foldid, times = d)
    } else if (length(foldid) == N * d) {
      foldid_glm <- foldid
    } else {
      warning("foldid must have length N or N*d; ignoring and letting glmnet choose folds.")
    }
  }
  
  # fit (note intercept = FALSE because we explicitly control intercept columns when desired)
  cvfit <- glmnet::cv.glmnet(
    x = X, y = y_stack,
    family = "gaussian",
    intercept = FALSE,
    standardize = standardize,
    penalty.factor = pf,
    foldid = foldid_glm,
    lambda = lambda
  )
  
  #lam_best <- cvfit$lambda.min  # or lambda.1se if you prefer
  if(lambest == "min"){
    
    lam_best <- cvfit$lambda.min
    
  } else  if(lambest == "1se") {
    
    lam_best <- cvfit$lambda.1se 
    
  } else {
  
    stop("lambest must be 'min' or '1se'")
    
  }
  
  
  # coef() is an S3 method; DON'T prefix with glmnet::  (that would error)
  b <- as.numeric(coef(cvfit, s = lam_best, exact = FALSE))
  # first entry is the (global) model intercept (always present, zero since intercept=FALSE)
  b <- b[-1L]
  
  if (include_intercept) {
    # expected length: d (per-outcome intercepts) + d*d (slopes)
    b0       <- b[seq_len(d)]
    bslopes  <- b[-seq_len(d)]
  } else {
    # no per-outcome intercepts; keep a zero column to preserve d x (d+1) shape
    b0       <- rep(0, d)
    bslopes  <- b
  }
  
  # arrange slopes into d x d where row = outcome, col = predictor
  slopes <- matrix(bslopes, nrow = d, ncol = d, byrow = TRUE)
  
  A_pre <- cbind(b0, slopes)  # d x (d+1)
  
  list(
    A_pre          = A_pre,
    lambda_best    = lam_best,
    cvfit          = cvfit,
    include_intercept = include_intercept
  )
}

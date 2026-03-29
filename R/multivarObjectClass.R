#comment
# Check multivar object.
check.multivar <- function(object){

    errors <- character()

    if(any(is.na(object@A))){
      msg <- c("multivar error: remove NA values before running constructModel.")
      errors <- c(errors,msg)
    }

    if(length(errors)==0) TRUE else errors

}

#' multivar object class
#'
#' An object class to be used with cv.multivar
#' 
#' @slot k Numeric. The number of subjects (or groupings) in the dataset.
#' @slot n Numeric Vector. Vector containing the number of timepoints for each dataset.
#' @slot d Numeric Vector. Vector containing the number of variables for each dataset.
#' @slot Ak List. A list (length = k) of lagged (T-lag-horizon) by d multivariate time series.
#' @slot bk List. A list (length = k) of (T-lag-horizon) by d multivariate time series.
#' @slot Ak_orig List. A list (length = k) of lagged (T-lag-horizon) by d unscaled multivariate time series.
#' @slot bk_orig List A list (length = k) of (T-lag-horizon) by d unscaled multivariate time series.
#' @slot Hk List. A list (length = k) of (horizon) by d multivariate time series.
#' @slot A Matrix. A matrix containing the lagged ((T-lag-horizon)k) by (d+dk) multivariate time series.
#' @slot b Matrix. A matrix containing the non-lagged ((T-lag-horizon)k) by (d) multivariate time series.
#' @slot H Matrix. A matrix containing the non-lagged (horizon k) by d multivariate time series.
#' @slot data_means List. A list (length = k) of column means for each group's data, used for intercept recovery when intercept=TRUE.
#' @slot lag Numeric. The VAR order. Currently only lag 1 is supported.
#' @slot horizon Numeric. Forecast horizon.
#' @slot t1 Numeric vector. Index of time series in which to start cross validation for individual k. 
#' @slot t2 Numeric vector. Index of time series in which to end cross validation for individual k.
#' @slot t1k Numeric vector. Index of time series in which to start cross validation for individual k.
#' @slot t2k Numeric vector. Index of time series in which to end cross validation for individual k.
#' @slot ntk Numeric. Number of usable timepoints (rows of b) per individual k
#' @slot ndk Numeric. Number of variables (cols of b) per individual k
#' @slot lambda1 Numeric vector. Regularization parameter grid.
#' @slot nlambda1 Numeric. Number of lambda1 values to search over. Default is 30.
#' @slot n_ratios_subgroup Numeric. Number of ratios_subgroup values to search over. Default is 30.
#' @slot gamma Numeric. Need definition here
#' @slot tol Numeric. Convergence tolerance.
#' @slot depth Numeric. Depth of grid construction. Default is 1000.
#' @slot window Numeric. Size of rolling window.
#' @slot standardize Logical. Default is true. Whether to standardize the individual data.
#' @slot weightest Character. How to estimate initial coefficients for adaptive weights. Default is "lasso". Other options include "ridge" and "ols". 
#' @slot canonical Logical. Default is false. If true, individual datasets are fit to a VAR(1) model.
#' @slot threshold Logical. Default is false. If true, and canonical is true, individual transition matrices are thresholded based on significance.
#' @slot lassotype Character. Default is "adaptive". Choices are "standard" or "adaptive" lasso.
#' @slot intercept Logical. Default is FALSE.
#' @slot W Matrix. Default is NULL. 
#' @slot ratios_unique Numeric vector. Penalty ratio for unique effects. Default is NULL.
#' @slot ratios_subgroup Numeric vector. Penalty ratio for subgroup effects. Default is NULL.
#' @slot ratios_unique_tvp Numeric vector. Default is NULL.
#' @slot ratios_common_tvp Numeric vector. Penalty ratio for common TVP effects. Default is NULL. 
#' @slot cv Character. Default is "blocked" for k-folds blocked cross-validation. rolling window cross-validation also available using "rolling".  If "blocked" is selected the nfolds argument should be specified.
#' @slot nfolds Numeric. The number of folds for use with "blocked" cross-validation.
#' @slot lamadapt Logical. Should the lambdas be calculated adaptively. Default is FALSE.
#' @slot subgroup_membership Numeric. Vector of subgroup assignments.
#' @slot subgroup Logical. Argument whether to run subgrouping algorithm.
#' @slot B Matrix. Default is NULL. 
#' @slot initcoefs List. A list of initial consistent estimates for the total, subgroup, unique and common effects.
#' @slot pendiag Logical. Logical indicating where autoregressive paramaters should be penalized. Default is true.
#' @slot tvp Logical.
#' @slot inittvpcoefs List. A list of initial tvp estimates.
#' @slot breaks List. A list of length K indicating structural breaks in the time series.
#' @slot lambda_choice Character. Which lambda to use for initial coefficient estimation ("lambda.min" or "lambda.1se").
#' @slot common_effects Logical. Whether to include common effects in TVP models. Only applies when tvp = TRUE.
#' @slot common_tvp_effects Logical. Whether to include common TVP effects (shared time-varying patterns) in TVP models. Only applies when tvp = TRUE.
#' @slot save_beta Logical. Whether to retain the full beta array in the cv.multivar result. Default is TRUE.
#' @slot ncores Numeric. Number of cores for parallel computation. Default is 1.
#' @slot spec List. Design matrix specification object created by \code{\link{build_matrix_spec}}. Single source of truth for column/row indices.
#' @slot eps Numeric. FISTA convergence tolerance. Default is 1e-3.
#' @slot warmstart Logical. Whether to use warm starts in the FISTA solver. Default is TRUE.
#' @slot stopping_crit Character. FISTA convergence criterion. One of \code{"absolute"}, \code{"relative"}, or \code{"objective"}.
#' @slot selection Character. Model selection criterion: \code{"cv"} (default) for cross-validated MSFE, or \code{"ebic"} for Extended BIC.
#' @slot ebic_gamma Numeric. EBIC tuning parameter, used when \code{selection = "ebic"}. Default is 0.5.
#' @slot weight_type Character. Adaptive weight function type: \code{"standard"} (default) uses \code{1/|coef|^gamma}, \code{"bounded"} uses \code{1/(1+|coef|/tau)^gamma}.
#' @details To construct an object of class multivar, use the function \code{\link{constructModel}}
#' @seealso \code{\link{constructModel}}
#' @export
setClass(
    Class = "multivar",
    representation(
        k = "numeric",
        n = "numeric",
        d = "numeric",
        Ak = "list",
        bk = "list",
        Ak_orig = "list",
        bk_orig = "list",
        Hk = "list",
        A = "matrix",
        b = "matrix",
        H = "matrix",
        data_means = "list",
        lag="numeric",
        horizon="numeric",
        t1 = "numeric",
        t2 = "numeric",
        t1k = "numeric",
        t2k = "numeric",
        ntk = "numeric",
        ndk = "numeric",
        lambda1 = "matrix",
        nlambda1 = "numeric",
        n_ratios_subgroup = "numeric",
        gamma = "numeric",
        tol = "numeric",
        depth = "numeric",
        window = "numeric",
        standardize = "logical",
        weightest = "character",
        canonical = "logical",
        threshold = "logical",
        lassotype = "character",
        intercept = "logical",
        W = "array",
        ratios_unique = "numeric",
        ratios_subgroup = "numeric",
        ratios_unique_tvp = "numeric",
        ratios_common_tvp = "numeric",
        cv = "character",
        nfolds = "numeric",
        lamadapt = "logical",
        subgroup_membership = "numeric",
        subgroup = "logical",
        B = "array",
        initcoefs = "list",
        pendiag = "logical",
        tvp = "logical",
        inittvpcoefs = "list",
        breaks = "list",
        lambda_choice = "character",
        common_effects = "logical",
        common_tvp_effects = "logical",
        save_beta = "logical",
        ncores = "numeric",
        spec = "list",
        eps = "numeric",
        warmstart = "logical",
        stopping_crit = "character",
        selection = "character",
        ebic_gamma = "numeric",
        weight_type = "character",
        maity_opts = "list"
        ),validity=check.multivar
    )




# show-default method to show an object when its name is printed in the console.
#' Default show method for an object of class multivar
#'
#' @param object \code{multivar} object created from \code{ConstructModel}
#' @return Displays the following information about the multivar object:
#' \itemize{
#' \item{To do.}
#' }
#' @seealso \code{\link{constructModel}} 
#' @name show.multivar
#' @aliases show,multivar-method
#' @docType methods
#' @rdname show-methods
#' @export
setMethod("show","multivar",function(object){
  cat("*** multivar model *** \n")
  cat("Number of groupings: ") ; cat(object@k, "\n")
  cat("Forecast horizon: ") ; cat(object@horizon, "\n")
  cat("Number of variables: \n") ; cat("    ");print(object@d)
  cat("Number of timepoints: \n") ; cat("    ");print(object@n)
})




#' Cross Validation for multivar
#' 
#' @usage cv.multivar(object)
#' @param object multivar object built using \code{ConstructModel}.
#' @details The main function of the multivar package. Performs cross validation to select penalty parameters over a training sample and evaluates them over a test set.
#' @return An object of class \code{multivar.results}.
#' @name cv.multivar
#' @aliases cv.multivar,multivar-method
#' @docType methods
#' @rdname cv.multivar-methods
#' @examples
#' 
#' # example 1 (run)
#' sim1  <- multivar_sim(
#'   k = 2,  # individuals
#'   d = 5,  # number of variables
#'   n = 20, # number of timepoints
#'   prop_fill_com = 0.1, # proportion of paths common
#'   prop_fill_ind = 0.05, # proportion of paths unique
#'   lb = 0.1,  # lower bound on coefficient magnitude
#'   ub = 0.5,  # upper bound on coefficient magnitude
#'   sigma = diag(5) # noise
#' )
#' 
#' model1 <- constructModel(data = sim1$data)
#' fit1 <- multivar::cv.multivar(model1)
#'
#'
#' @export
setGeneric(name = "cv.multivar",def=function(object){standardGeneric("cv.multivar")})
setMethod(f = "cv.multivar", signature = "multivar",definition = function(object){

  object@initcoefs <- estimate_initial_coefs(
    object@Ak,
    object@bk,
    object@d,
    object@k,
    object@lassotype,
    object@weightest,
    object@subgroup_membership,
    object@subgroup,
    object@nlambda1,
    object@tvp,
    object@breaks,
    object@intercept,
    object@nfolds,
    object@lambda_choice,
    object@common_effects,
    object@common_tvp_effects,
    object@maity_opts
  )

  object@W <- est_base_weight_mat(
    object@W,
    object@Ak,
    object@initcoefs,
    object@ratios_unique,
    object@d,
    object@k,
    object@lassotype,
    object@weightest,
    object@subgroup_membership,
    object@subgroup,
    object@ratios_subgroup,
    object@pendiag,
    object@tvp,
    object@ratios_unique_tvp,
    object@ratios_common_tvp,
    object@intercept,
    object@common_effects,
    object@common_tvp_effects,
    object@spec,
    object@weight_type
  )

  # Precompute transposes once (used for lambda_grid and cv_multivar)
  tA <- t(as.matrix(object@A))
  tb <- t(as.matrix(object@b))

  # Only construct lambda grid if user didn't provide fixed lambda values
  # User-provided lambdas have non-zero values; default is all zeros
  if (all(object@lambda1 == 0)) {
    object@lambda1 <- lambda_grid(
      depth     = object@depth,
      nlam      = object@nlambda1,
      Y         = tb,
      Z         = tA,
      W         = object@W,
      tol       = object@tol,
      intercept = object@intercept,
      lamadapt  = object@lamadapt,
      k         = object@k
    )
  }

  stopping_crit_int <- match(object@stopping_crit, c("absolute", "relative", "objective")) - 1L

  if (object@selection == "ebic") {
    # EBIC path: fit all scenarios on full data (no CV folds)
    beta <- wlasso(object@B, tA, tb, object@W, object@k, object@d,
                   object@lambda1, object@eps, object@intercept,
                   warmstart = object@warmstart, stopping_crit = stopping_crit_int)
    fit <- list(beta, NULL)  # No MSFE matrix

    # Select best scenario by EBIC
    ebic_vals <- compute_ebic(beta, tA, tb, object@d, object@ebic_gamma, object@spec)
    best_idx <- which.min(ebic_vals)

  } else {
    # CV path: existing cross-validation code
    fit <- cv_multivar(
      object@B,
      tA,
      tb,
      object@W,
      object@Ak,
      object@bk,
      object@k,
      object@d,
      object@lambda1,
      object@t1,
      object@t2,
      eps = object@eps,
      object@intercept,
      object@cv,
      object@nfolds,
      object@tvp,
      object@breaks,
      object@spec,
      object@ncores,
      warmstart = object@warmstart,
      stopping_crit = stopping_crit_int
    )

    hyp <- extract_multivar_hyperparams(
      object,
      fit
    )

    # Check if selected hyperparameters are at grid boundaries
    best_idx <- which.min(hyp$MSFE)
    best_lambda_idx <- hyp$lambda1_index[best_idx]
    best_ratio_idx <- hyp$ratio_index[best_idx]
    n_lambda <- nrow(object@lambda1)
    n_ratios <- ncol(object@lambda1)

    # Check lambda1 boundaries
    if (best_lambda_idx == 1) {
      warning("lambda1 selected at upper boundary (maximum regularization). ",
              "Consider increasing depth or checking if model is overpenalized.")
    } else if (best_lambda_idx == n_lambda) {
      warning("lambda1 selected at lower boundary (minimum regularization). ",
              "Consider increasing depth for wider lambda range.")
    }

    # Check ratio boundaries (only warn if more than 1 ratio value)
    if (n_ratios > 1) {
      if (best_ratio_idx == 1) {
        warning("ratio selected at upper boundary. ",
                "Consider expanding the ratios_unique grid.")
      } else if (best_ratio_idx == n_ratios) {
        warning("ratio selected at lower boundary. ",
                "Consider expanding the ratios_unique grid.")
      }
    }
  }

  mats <- breakup_transition(
    B = fit[[1]][,,best_idx],
    spec = object@spec,
    Ak = object@Ak,
    breaks = object@breaks
  )

  # Recover intercepts if intercept=TRUE
  # Intercepts are computed post-hoc using: c = mean(b) - Phi * mean(A)
  if (object@intercept && length(object@data_means) > 0) {
    mats$intercepts <- recover_intercepts(
      mats = mats,
      data_means = object@data_means,
      k = object@k,
      d = object@d,
      subgroup = object@subgroup,
      subgroup_membership = object@subgroup_membership,
      tvp = object@tvp,
      breaks = object@breaks
    )
  }

  if (object@selection == "ebic") {
    results <- list(
      mats = mats,
      beta = if (object@save_beta) fit[[1]] else NULL,
      MSFE = NULL,
      obj  = object,
      hyperparams = NULL,
      selection = "ebic",
      ebic = list(
        values   = ebic_vals,
        best_idx = best_idx,
        gamma    = object@ebic_gamma
      )
    )
  } else {
    results <- list(
      mats = mats,
      beta = if (object@save_beta) fit[[1]] else NULL,
      MSFE = fit[[2]],
      obj  = object,
      hyperparams = hyp,
      selection = "cv"
    )
  }

  # Add S3 class for method dispatch
  class(results) <- c("multivar_fit", "list")

  #results <- new("multivar.results",object)
  return(results)
})

#' Canonical VAR Fitting Function for multivar
#' 
#' @usage canonical.multivar(object)
#' @param object multivar object built using \code{ConstructModel}.
#' @details A function to fit a canonical VAR model to each individual dataset. 
#' @return A list of results.
#' @seealso \code{\link{constructModel}}, 
#' @name canonical.multivar
#' @aliases canonical.multivar,multivar-method
#' @docType methods
#' @rdname canonical.multivar-methods
#' @examples
#' 
#' # example 1 (run)
#' sim1  <- multivar_sim(
#'   k = 2,  # individuals
#'   d = 5,  # number of variables
#'   n = 20, # number of timepoints
#'   prop_fill_com = 0.1, # proportion of paths common
#'   prop_fill_ind = 0.05, # proportion of paths unique
#'   lb = 0.1,  # lower bound on coefficient magnitude
#'   ub = 0.5,  # upper bound on coefficient magnitude
#'   sigma = diag(5) # noise
#' )
#' 
#' model1 <- constructModel(data = sim1$data, weightest = "ols")
#' fit1 <- canonical.multivar(model1)
#'
#' @export
setGeneric(name = "canonical.multivar",def=function(object){standardGeneric("canonical.multivar")})
setMethod(f = "canonical.multivar", signature = "multivar",definition = function(object){

  can_var <- lapply(seq_along(object@bk_orig), function(i){
    # Use original (unscaled) data for canonical VAR fitting
    # Reconstruct the time series by prepending one lag observation
    df <- as.matrix(rbind(
      object@Ak_orig[[i]][1,],
      object@bk_orig[[i]]
     ))
     # Use type="const" if intercept=TRUE, otherwise type="none"
     var_type <- if(object@intercept) "const" else "none"
     fit_canonical_var(df, p = 1, type = var_type)
  })


  res <- list(
    common = NULL,
    unique = NULL,
    total  = lapply(can_var,"[[","transition_mat"),
    total_sigonly  = lapply(can_var,"[[","transition_mat_sigonly")
  )

  # Handle intercepts
  if (object@intercept) {
    # Extract intercepts directly from canonical VAR fit (vars with type="const")
    intercepts_total <- lapply(can_var, "[[", "intercepts")

    # For k=1, no decomposition
    if (object@k == 1) {
      res$intercepts <- list(
        intercepts_total = intercepts_total,
        intercept_common = NULL,
        intercepts_unique = NULL
      )
    } else {
      # For k>1, decompose into common and group-specific
      intercept_common <- Reduce("+", intercepts_total) / object@k
      intercepts_unique <- lapply(intercepts_total, function(c_total) {
        c_total - intercept_common
      })

      res$intercepts <- list(
        intercepts_total = intercepts_total,
        intercept_common = intercept_common,
        intercepts_unique = intercepts_unique
      )
    }
  } else {
    res$intercepts <- NULL
  }

  return(list(mats = res))
})














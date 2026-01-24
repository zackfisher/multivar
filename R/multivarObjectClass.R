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
#' @slot Hk List. A list (length = k) of (horizon) by d multivariate time series.
#' @slot A  Matrix. A matrix containing the lagged ((T-lag-horizon)k) by (d+dk) multivariate time series.
#' @slot b  Matrix. A matrix containing the non-lagged ((T-lag-horizon)k) by (d) multivariate time series.
#' @slot H  Matrix. A matrix containing the non-lagged (horizon k) by d multivariate time series.
#' @slot lag Numeric. The VAR order. Currently only lag 1 is supported.
#' @slot horizon Numeric. Forecast horizon.
#' @slot t1 Numeric vector. Index of time series in which to start cross validation for individual k. 
#' @slot t2 Numeric vector. Index of time series in which to end cross validation for individual k.
#' @slot lambda1 Numeric vector. Regularization parameter 1.
#' @slot lambda2 Numeric vector. Regularization parameter 2.
#' @slot tau Numeric vector. Regularization parameter for subgroup effects.
#' @slot nlambda1 Numeric. Number of lambda1 values to search over. Default is 30.
#' @slot nlambda2 Numeric. Number of lambda2 values to search over. Default is 30.
#' @slot ntau Numeric. Number of tau values to search over. Default is 30.
#' @slot tol Numeric. Convergence tolerance.
#' @slot depth Numeric. Depth of grid construction. Default is 1000.
#' @slot window Numeric. Size of rolling window.
#' @slot standardize Logical. Default is true. Whether to standardize the individual data.
#' @slot weightest Character. How to estimate the first-stage weights. Default is "lasso". Other options include "ridge" and "ols". 
#' @slot canonical Logical. Default is false. If true, individual datasets are fit to a VAR(1) model.
#' @slot threshold Logical. Default is false. If true, and canonical is true, individual transition matrices are thresholded based on significance.
#' @slot lassotype Character. Default is "adaptive". Choices are "standard" or "adaptive" lasso.
#' @slot intercept Logical. Default is FALSE.
#' @slot pen_common_intercept Logical. Default is FALSE. Whether to penalize the common intercept.
#' @slot pen_unique_intercept Logical. Default is TRUE. Whether to penalize individual-specific intercept deviations.
#' @slot W Matrix. Default is NULL. 
#' @slot ratios Numeric vector. Default is NULL. 
#' @slot ratiostau Numeric vector. Default is NULL. 
#' @slot ratiosalpha Numeric vector. Default is NULL. 
#' @slot cv Character. Default is "blocked" for k-folds blocked cross-validation. rolling window cross-validation also available using "rolling".  If "blocked" is selected the nfolds argument should be specified.
#' @slot nfolds Numeric. The number of folds for use with "blocked" cross-validation.
#' @slot thresh Numeric. Post-estimation threshold for setting the individual-level coefficients to zero if their absolute value is smaller than the value provided. Default is zero.
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
#' @details To construct an object of class multivar, use the function \code{\link{constructModel}}
#' @seealso \code{\link{constructModel}}
#' @export
setClass(
    Class="multivar",
    representation(
        k = "numeric",
        n  = "numeric",
        d  = "numeric",
        Ak = "list",
        bk = "list",
        Hk = "list",
        A  = "matrix",
        b  = "matrix",
        H  = "matrix",
        lag="numeric",
        horizon="numeric",
        t1 = "numeric",
        t2 = "numeric",
        t1k = "numeric",
        t2k = "numeric",
        ntk = "numeric",
        ndk = "numeric",
        lambda1="matrix",
        lambda2="matrix",
        tau="matrix",
        nlambda1="numeric",
        nlambda2="numeric",
        ntau ="numeric",
        gamma = "numeric",
        tol="numeric",
        depth="numeric",
        window="numeric",
        standardize = "logical",
        weightest = "character",
        canonical = "logical",
        threshold = "logical",
        lassotype = "character",
        intercept = "logical",
        pen_common_intercept = "logical",
        pen_unique_intercept = "logical",
        W = "array",
        ratios = "numeric",
        ratiostau = "numeric",
        ratiosalpha = "numeric",
        cv = "character",
        nfolds = "numeric",
        thresh = "numeric",
        lamadapt = "logical",
        subgroup_membership = "numeric",
        subgroup = "logical",
        B = "array",
        initcoefs = "list",
        pendiag = "logical",
        tvp = "logical",
        inittvpcoefs = "list",
        breaks = "list",
        lambda_choice = "character"
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
    object@nlambda2,
    object@tvp,
    object@breaks,
    object@intercept,
    object@nfolds,
    object@lambda_choice
  )
  
  object@W <- est_base_weight_mat(
    object@W,
    object@Ak,
    object@initcoefs,
    object@ratios,
    object@d,
    object@k,
    object@lassotype,
    object@weightest,
    object@subgroup_membership,
    object@subgroup,
    object@ratiostau,
    object@pendiag,
    object@tvp,
    object@ratiosalpha,
    object@intercept,
    object@pen_common_intercept,
    object@pen_unique_intercept
  )

  # object@W <- est_base_weight_mat(
  #   object@W,
  #   object@Ak,
  #   object@bk,
  #   object@ratios, 
  #   object@d, 
  #   object@k, 
  #   object@lassotype, 
  #   object@weightest,
  #   object@subgroup,
  #   object@subgroupflag,
  #   object@ratiostau,
  #   object@nlambda1,
  #   object@nlambda2
  # )
  
  object@lambda1 <- lambda_grid(
    #object@B, 
    object@depth, 
    object@nlambda1, 
    t(as.matrix(object@b)), 
    t(as.matrix(object@A)), 
    object@W, 
    object@k,
    object@tol,
    object@intercept,
    object@lamadapt
  )

  fit <- cv_multivar(
    object@B, 
    t(as.matrix(object@A)), 
    t(as.matrix(object@b)), 
    object@W, 
    object@Ak,
    object@bk,
    object@k, 
    object@d, 
    object@lambda1, 
    object@t1, 
    object@t2, 
    eps = 1e-3,
    object@intercept,
    object@cv,
    object@nfolds
  )
  

  hyp <- extract_multivar_hyperparams(
    object, 
    fit
  )


  mats <- breakup_transition(
    fit[[1]][,,which.min(colMeans(fit[[2]]))], 
    object@Ak, 
    object@ndk, 
    object@intercept,
    object@thresh,
    object@subgroup_membership,
    object@subgroup,
    object@tvp,
    object@ntk,
    object@breaks
  )
  
  results <- list(
    mats = mats,
    beta = fit[[1]],
    MSFE = fit[[2]],
    obj  = object,
    hyperparams = hyp
  )
  
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

  can_var <- lapply(seq_along(object@bk), function(i){
    df <- as.matrix(rbind(
      object@Ak[[i]][1,],
      object@bk[[i]]
     ))
     fit_canonical_var(df, p = 1, type = "none")
  })


  res <- list(
    common = NULL,
    unique = NULL,
    total  = lapply(can_var,"[[","transition_mat"),
    total_sigonly  = lapply(can_var,"[[","transition_mat_sigonly")
  )

  return(list(mats = res))
})














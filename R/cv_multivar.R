#' Cross-Validation Dispatcher for multivar
#'
#' Routes to the appropriate CV method (blocked or rolling).
#'
#' @param B Initial coefficient array
#' @param Z Design matrix (transposed)
#' @param Y Response matrix (transposed)
#' @param W Weight array
#' @param Ak List of design matrices per subject
#' @param bk List of response matrices per subject
#' @param k Number of subjects
#' @param d Number of variables
#' @param lambda1 Lambda grid
#' @param t1 Start indices
#' @param t2 End indices
#' @param eps Convergence tolerance
#' @param intercept Whether model includes intercepts
#' @param cv CV method: "blocked" or "rolling"
#' @param nfolds Number of CV folds
#' @param tvp Whether time-varying parameters are used
#' @param breaks List of period indices per subject
#' @param spec Optional matrix_spec object for row/column indices
#'
#' @return List with beta coefficients and MSFE matrix
#' @export
cv_multivar <- function(B, Z, Y, W, Ak, bk, k, d, lambda1, t1, t2, eps,
                        intercept = FALSE, cv, nfolds, tvp = FALSE,
                        breaks = NULL, spec = NULL) {

  if (cv == "rolling") {
    # TODO: Add spec support to cv_rolling when needed
    res <- cv_rolling(B, Z, Y, W, Ak, k, d, lambda1, t1, t2, eps, intercept, cv, nfolds)

  } else if (cv == "blocked") {
    res <- cv_blocked(B, Z, Y, W, Ak, k, d, lambda1, t1, t2, eps,
                      intercept, cv, nfolds, tvp, breaks, spec)

  } else {
    stop(paste0("multivar ERROR: ", cv, " is not a supported cross-validation method."))
  }

  return(res)
}
#' Extract hyperparameter grid search results
#'
#' Creates a tidy dataframe with all tested hyperparameter combinations and their
#' cross-validation error (MSFE). Each row represents one (lambda1, ratio) combination.
#'
#' @param object multivar model object used to call cv_multivar(...)
#' @param fit    result returned by cv_multivar(...); MSFE is assumed to be fit[[2]]
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{ratio_index}{Index of ratio (1 to n_ratios)}
#'     \item{lambda1_index}{Index of lambda1 (1 to n_lambda)}
#'     \item{slice_index}{Flattened index using C++ ordering}
#'     \item{lambda1_value}{Actual lambda1 penalty value}
#'     \item{ratio_value}{Actual ratio value}
#'     \item{lambda2_value}{Computed lambda2 = lambda1 * ratio}
#'     \item{MSFE}{Mean squared forecast error averaged across CV folds}
#'   }
#'
#' @details
#' The function averages MSFE across cross-validation folds and returns one row
#' per hyperparameter combination. The dataframe can be used for:
#' \itemize{
#'   \item Finding best hyperparameters: \code{df[which.min(df$MSFE), ]}
#'   \item Plotting MSFE surface
#'   \item Filtering/analyzing hyperparameter performance
#' }
#'
#' @export
extract_multivar_hyperparams <- function(object, fit) {

  # Number of lambda1s and ratios
  n_lambda <- nrow(object@lambda1)
  n_ratios <- ncol(object@lambda1)

  # Create all combinations of lambda1 and ratio indices
  param_map <- expand.grid(
    ratio_index   = seq_len(n_ratios),
    lambda1_index = seq_len(n_lambda)
  )

  # Compute slice index consistent with C++ ordering: i*n_r + j
  param_map$slice_index <- (param_map$lambda1_index - 1) * n_ratios + (param_map$ratio_index - 1)

  # Assign the actual lambda1 and ratio values
  param_map$lambda1_value <- as.vector(t(object@lambda1))
  param_map$ratio_value   <- object@ratios_unique[param_map$ratio_index]

  # Compute lambda2 = lambda1 * ratio
  param_map$lambda2_value <- param_map$lambda1_value * param_map$ratio_value

  # Extract MSFE and average across folds
  # Handle both 2D (nfolds x n_combos) and 3D (d x ? x n_combos) arrays
  msfe_raw <- fit[[2]]

  if (length(dim(msfe_raw)) == 2) {
    # 2D case: (nfolds x n_combos)
    MSFE_vec <- colMeans(msfe_raw)
  } else if (length(dim(msfe_raw)) == 3) {
    # 3D case: average across first two dimensions
    MSFE_vec <- apply(msfe_raw, 3, mean)
  } else {
    stop("extract_multivar_hyperparams: unexpected MSFE array dimensions: ",
         paste(dim(msfe_raw), collapse = " x "))
  }

  # Verify correct length
  expected_length <- n_lambda * n_ratios
  if (length(MSFE_vec) != expected_length) {
    stop("extract_multivar_hyperparams: MSFE length mismatch. ",
         "Expected ", expected_length, ", got ", length(MSFE_vec))
  }

  param_map$MSFE <- MSFE_vec

  return(param_map)
}

#' Compute EBIC for an array of fitted coefficient matrices
#'
#' Core computation of the Extended Bayesian Information Criterion (EBIC)
#' for each scenario in a beta array. Used internally by \code{cv.multivar}
#' when \code{selection = "ebic"} and by \code{select_by_ebic} for post-hoc
#' reselection.
#'
#' @param beta_array 3D array (d x p x n_scenarios) of fitted coefficients.
#' @param Z Design matrix (transposed), p x n.
#' @param Y Response matrix (transposed), d x n.
#' @param d Number of response variables.
#' @param gamma EBIC tuning parameter. \code{gamma = 0} gives standard BIC;
#'   \code{gamma = 0.5} (default) is a moderate EBIC; \code{gamma = 1} is
#'   the most conservative (strongest sparsity preference).
#'
#' @return Numeric vector of EBIC values, one per scenario.
#'
#' @details
#' The EBIC is computed per-equation and summed:
#' \deqn{\text{EBIC}_i = \sum_{j=1}^{d} n \log(\text{RSS}_{ij} / n) + k_i \log(n) + 2 \gamma \, k_i \log(p)}
#' where \eqn{\text{RSS}_{ij}} is the residual sum of squares for equation
#' \eqn{j} under scenario \eqn{i}, \eqn{k_i} is the total number of
#' nonzero coefficients across all equations, and \eqn{p} is the number
#' of predictors per equation.
#'
#' @export
compute_ebic <- function(beta_array, Z, Y, d, gamma = 0.5) {
  n <- ncol(Z)
  p <- nrow(Z)
  n_scenarios <- dim(beta_array)[3]
  ebic_vals <- numeric(n_scenarios)

  for (i in seq_len(n_scenarios)) {
    B_i <- beta_array[, , i]
    resid <- Y - B_i %*% Z
    rss_per_eq <- pmax(rowSums(resid^2), .Machine$double.eps)
    k_i <- sum(abs(B_i) > 0)
    ebic_vals[i] <- sum(n * log(rss_per_eq / n)) + k_i * log(n) + 2 * gamma * k_i * log(p)
  }

  ebic_vals
}


#' Reselect best model using (E)BIC instead of MSFE
#'
#' Takes the output of \code{cv.multivar()} and reselects the best penalty
#' scenario using the Extended Bayesian Information Criterion (EBIC) or
#' standard BIC. This can improve structure recovery (fewer false positives)
#' compared to prediction-based cross-validation, especially when n >> p.
#'
#' @param fit Output of \code{cv.multivar()}.
#' @param gamma EBIC tuning parameter. \code{gamma = 0} gives standard BIC;
#'   \code{gamma = 0.5} (default) is a moderate EBIC; \code{gamma = 1} is
#'   the most conservative (strongest sparsity preference).
#'
#' @return A modified copy of \code{fit} with:
#'   \item{mats}{Transition matrices from the EBIC-selected scenario}
#'   \item{ebic}{Named list: \code{values} (EBIC for all scenarios),
#'         \code{best_idx} (selected scenario index), \code{gamma} (tuning value)}
#'   The original MSFE results are preserved.
#'
#' @export
select_by_ebic <- function(fit, gamma = 0.5) {

  if (is.null(fit$beta)) {
    stop("select_by_ebic requires the full beta array. ",
         "Re-run cv.multivar() with save_beta = TRUE (the default).")
  }

  obj <- fit$obj

  # Reconstruct Z and Y (transposed, as used in cv_blocked)
  Z <- t(as.matrix(obj@A))
  Y <- t(as.matrix(obj@b))

  ebic_vals <- compute_ebic(fit$beta, Z, Y, obj@d, gamma)
  best_idx <- which.min(ebic_vals)

  # Rebuild mats from the EBIC-selected scenario
  mats <- breakup_transition(
    B = fit$beta[, , best_idx],
    spec = obj@spec,
    Ak = obj@Ak,
    breaks = obj@breaks
  )

  # Recover intercepts if needed
  if (obj@intercept && length(obj@data_means) > 0) {
    mats$intercepts <- recover_intercepts(
      mats = mats,
      data_means = obj@data_means,
      k = obj@k,
      d = obj@d,
      subgroup = obj@subgroup,
      subgroup_membership = obj@subgroup_membership,
      tvp = obj@tvp,
      breaks = obj@breaks
    )
  }

  # Return modified fit
  fit$mats <- mats
  fit$ebic <- list(
    values   = ebic_vals,
    best_idx = best_idx,
    gamma    = gamma
  )

  fit
}

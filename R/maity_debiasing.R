#' Maity et al. (2022) MrLasso-style initial estimation
#'
#' Provides debiased, robustly-aggregated, analytically-thresholded initial
#' estimates for adaptive LASSO weights. Replaces the default median-based
#' decomposition when \code{weightest = "maity"} in \code{\link{constructModel}}.
#'
#' @references
#' Maity, S., Sun, Y., & Banerjee, M. (2022). Meta-analysis of heterogeneous
#' data: integrative sparse regression in high-dimensions. \emph{JMLR}, 23, 1-50.
#'
#' @name maity_debiasing
NULL


#' Estimate initial coefficients using Maity et al. (2022) MrLasso approach
#'
#' Five-step pipeline:
#' \enumerate{
#'   \item Per-subject LASSO via \code{\link{estimate_raw_effects}}
#'   \item Debias each subject via nodewise regression
#'   \item Aggregate common effects via re-descending loss
#'   \item Soft-threshold common at rate sqrt(log(d^2) / (k * min(n_k)))
#'   \item Unique = debiased - common_dense, thresholded at sqrt(log(d^2) / n_k)
#' }
#'
#' @param Ak List of design matrices per subject (each n_k x d)
#' @param bk List of response matrices per subject (each n_k x d)
#' @param d Number of variables
#' @param k Number of subjects
#' @param nfolds Number of CV folds for LASSO fits
#' @param lambda_choice Which lambda from cv.glmnet ("lambda.min" or "lambda.1se")
#' @param eta Numeric. Scale parameter for re-descending loss. Default 2.0.
#' @param thresh_const Numeric. Multiplicative constant for threshold rates. Default 1.0.
#'
#' @return List with same structure as \code{\link{decompose_effects}}
#' @keywords internal
estimate_maity_effects <- function(Ak, bk, d, k, nfolds, lambda_choice,
                                    eta = 2.0, thresh_const = 1.0) {

  # Step 1: Per-subject LASSO estimates
  raw <- estimate_raw_effects(Ak, bk, weightest = "lasso", tvp = FALSE,
                               breaks = NULL, nfolds = nfolds,
                               lambda_choice = lambda_choice)

  # Step 2: Debias each subject's estimates
  debiased <- debias_var1_subjects(Ak, bk, raw$total_effects, nfolds = nfolds)

  # Step 3: Robust aggregation via re-descending loss
  common_dense <- aggregate_common_redescending(debiased, d, k, eta)

  # Step 4: Threshold common effects
  nk <- vapply(Ak, nrow, numeric(1))
  t_global <- thresh_const * sqrt(log(d^2) / (k * min(nk)))
  common_effects <- soft_threshold_matrix(common_dense, t_global)

  # Step 5: Unique deviations (from dense common), thresholded
  unique_effects <- vector("list", k)
  total_effects <- vector("list", k)
  for (g in seq_len(k)) {
    t_local <- thresh_const * sqrt(log(d^2) / nk[g])
    delta_g <- debiased[[g]] - common_dense
    unique_effects[[g]] <- soft_threshold_matrix(delta_g, t_local)
    total_effects[[g]] <- common_effects + unique_effects[[g]]
  }

  list(
    common_effects = common_effects,
    subgroup_effects = NULL,
    unique_effects = unique_effects,
    tvp_effects = NULL,
    common_tvp_effects = NULL,
    total_effects = total_effects
  )
}


#' Debias VAR(1) estimates for all subjects via nodewise regression
#'
#' For each subject, builds a nodewise LASSO precision estimate and applies
#' the debiasing correction equation-by-equation.
#'
#' @param Ak List of design matrices (each n_k x d)
#' @param bk List of response matrices (each n_k x d)
#' @param coef_list List of d x d LASSO coefficient matrices
#' @param nfolds Number of CV folds for nodewise regressions
#'
#' @return List of d x d debiased coefficient matrices
#' @keywords internal
debias_var1_subjects <- function(Ak, bk, coef_list, nfolds = 5) {
  k <- length(coef_list)
  debiased <- vector("list", k)

  for (g in seq_len(k)) {
    X_g <- Ak[[g]]
    Y_g <- bk[[g]]
    Phi_hat <- coef_list[[g]]
    debiased[[g]] <- debias_var1_single(X_g, Y_g, Phi_hat, nfolds)
  }

  debiased
}


#' Debias a single subject's VAR(1) coefficient matrix
#'
#' @param X Design matrix (n x d)
#' @param Y Response matrix (n x d)
#' @param Phi_hat LASSO estimate (d x d, row i = equation i)
#' @param nfolds Number of CV folds
#'
#' @return Debiased d x d coefficient matrix
#' @keywords internal
debias_var1_single <- function(X, Y, Phi_hat, nfolds = 5) {
  n <- nrow(X)
  d <- ncol(X)

  # Build nodewise LASSO precision estimate (d x d)
  Theta_hat <- build_theta_nodewise(X, nfolds)

  # Residuals under original LASSO: R = Y - X %*% t(Phi_hat)
  R <- Y - X %*% t(Phi_hat)

  # Debias each equation (column of the transposed system)
  Phi_tilde <- matrix(0, d, d)
  for (i in seq_len(d)) {
    beta_hat_i <- Phi_hat[i, ]
    r_i <- R[, i]
    score_i <- crossprod(X, r_i) / n
    Phi_tilde[i, ] <- beta_hat_i + as.vector(Theta_hat %*% score_i)
  }

  Phi_tilde
}


#' Build nodewise LASSO precision estimate
#'
#' For each column j of X, regress X_j on X_{-j} via LASSO.
#' Construct the precision matrix Omega from the nodewise coefficients
#' and residual variances, then symmetrize.
#'
#' @param X Design matrix (n x d)
#' @param nfolds Number of CV folds
#'
#' @return d x d precision estimate (symmetrized)
#' @keywords internal
build_theta_nodewise <- function(X, nfolds = 5) {
  n <- nrow(X)
  d <- ncol(X)
  Omega <- matrix(0, d, d)

  for (j in seq_len(d)) {
    y_j <- X[, j]
    X_minus_j <- X[, -j, drop = FALSE]

    cvfit <- glmnet::cv.glmnet(
      x = X_minus_j, y = y_j,
      intercept = FALSE, standardize = FALSE,
      alpha = 1, nfolds = nfolds
    )

    beta_j <- as.numeric(glmnet::coef.glmnet(
      cvfit$glmnet.fit, s = cvfit$lambda.min
    ))[-1]

    res_j <- y_j - X_minus_j %*% beta_j
    tau2_j <- sum(res_j^2) / n

    Omega[j, j] <- 1 / tau2_j
    Omega[-j, j] <- -beta_j / tau2_j
  }

  # Symmetrize
  0.5 * (Omega + t(Omega))
}


#' Coordinate-wise re-descending loss aggregation
#'
#' Apply \code{\link{redescending_aggregate}} element-wise across k
#' d x d matrices to produce a robust common estimate.
#'
#' @param debiased_list List of k d x d debiased matrices
#' @param d Number of variables
#' @param k Number of subjects
#' @param eta Scale parameter for re-descending loss
#'
#' @return d x d common effects matrix (dense, before thresholding)
#' @keywords internal
aggregate_common_redescending <- function(debiased_list, d, k, eta) {
  common <- matrix(0, d, d)

  for (i in seq_len(d)) {
    for (j in seq_len(d)) {
      values <- vapply(debiased_list, function(m) m[i, j], numeric(1))
      common[i, j] <- redescending_aggregate(values, eta)
    }
  }

  common
}


#' Scalar re-descending loss minimizer
#'
#' Finds mu that minimizes sum_k min((x_k - mu)^2, eta^2) via iterative
#' trimmed mean. Adapted from Maity et al. (2022) MrLasso.
#'
#' @param x Numeric vector of values across subjects
#' @param eta Scale parameter for re-descending loss
#' @param tol Convergence tolerance
#' @param max_iter Maximum iterations
#'
#' @return Scalar robust aggregate
#' @keywords internal
redescending_aggregate <- function(x, eta, tol = 1e-9, max_iter = 1000) {
  mu <- mean(x)

  for (iter in seq_len(max_iter)) {
    # Keep inliers: subjects within eta + |mu| of current estimate
    diffs <- abs(x - mu)
    inliers <- x[diffs <= eta + abs(mu)]

    if (length(inliers) == 0) {
      return(mu)
    }

    mu_new <- mean(inliers)
    if (abs(mu_new - mu) < tol) break
    mu <- mu_new
  }

  mu
}


#' Element-wise soft thresholding of a matrix
#'
#' @param mat Matrix to threshold
#' @param threshold Threshold value
#'
#' @return Thresholded matrix
#' @keywords internal
soft_threshold_matrix <- function(mat, threshold) {
  sign(mat) * pmax(abs(mat) - threshold, 0)
}

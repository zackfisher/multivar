#' Estimate initial coefficients for adaptive weights
#'
#' @description
#' `estimate_initial_coefs` returns consistent initial estimates to be used
#' for calculating adaptive weights in adaptive LASSO regularization.
#'
#' @param Ak List. A list (length = k) of lagged (T-lag-horizon) by d multivariate time series matrices.
#' @param bk List. A list (length = k) of (T-lag-horizon) by d multivariate time series matrices.
#' @param d Numeric vector. The number of variables for each dataset.
#' @param k Numeric. The number of subjects.
#' @param lassotype Character. Type of LASSO penalty: "standard" or "adaptive".
#' @param weightest Character. How to estimate initial coefficients: "lasso", "ridge", or "ols".
#' @param subgroup_membership Numeric vector. Vector of subgroup assignments for each subject.
#' @param subgroup Logical. Whether to run subgrouping algorithm.
#' @param nlambda1 Numeric. Not used (kept for API compatibility).
#' @param tvp Logical. Whether to estimate time-varying parameters.
#' @param breaks List. A list of length k indicating structural breaks in the time series.
#' @param intercept Logical. Whether to include intercepts in the model.
#' @param nfolds Numeric. The number of folds for cross-validation.
#' @param lambda_choice Character. Which lambda to use from cv.glmnet: "lambda.min" or "lambda.1se".
#' @param common_effects Logical. Whether to include common effects in TVP models.
#' @param common_tvp_effects Logical. Whether to include common TVP effects (k>1 only).
#'
#' @return A list containing:
#'   \itemize{
#'     \item{common_effects: d x d matrix of common effects}
#'     \item{subgroup_effects: list of d x d matrices per subgroup (if subgroup=TRUE)}
#'     \item{unique_effects: list of d x d matrices per subject}
#'     \item{tvp_effects: list of lists of d x d matrices per subject per period (if tvp=TRUE)}
#'     \item{common_tvp_effects: list of d x d matrices per period (if tvp=TRUE and k>1)}
#'     \item{total_effects: list of d x d matrices per subject}
#'   }
#'
#' @export
estimate_initial_coefs <- function(
    Ak,
    bk,
    d,
    k,
    lassotype,
    weightest,
    subgroup_membership,
    subgroup,
    nlambda1,
    tvp,
    breaks,
    intercept,
    nfolds,
    lambda_choice = "lambda.min",
    common_effects = TRUE,
    common_tvp_effects = TRUE,
    maity_opts = list()) {


  # Standard lasso: return all-ones matrices (no adaptive weighting)
  if (lassotype == "standard") {
    return(make_unit_weights(d, k, subgroup, subgroup_membership, tvp, Ak, breaks))
  }

  # Maity et al. (2022) MrLasso approach: debias + re-descending loss + threshold
  if (weightest == "maity") {
    eta <- if (!is.null(maity_opts$eta)) maity_opts$eta else 2.0
    thresh_const <- if (!is.null(maity_opts$thresh_const)) maity_opts$thresh_const else 1.0
    return(estimate_maity_effects(
      Ak = Ak, bk = bk, d = d, k = k,
      nfolds = nfolds, lambda_choice = lambda_choice,
      eta = eta, thresh_const = thresh_const
    ))
  }

  # Step 1: Estimate raw effects (per-subject, and per-period if TVP)
  raw_estimates <- estimate_raw_effects(
    Ak, bk, weightest, tvp, breaks, nfolds, lambda_choice
  )

  # Step 2: Decompose into common/unique/subgroup/tvp structure
  decompose_effects(
    raw_estimates, k, d, subgroup, subgroup_membership,
    tvp, common_effects, common_tvp_effects
  )
}


#' Create unit weight matrices for standard (non-adaptive) LASSO
#'
#' @keywords internal
make_unit_weights <- function(d, k, subgroup, subgroup_membership, tvp, Ak, breaks) {

  if (tvp) {
    list(
      common_effects = matrix(1, d, d),
      subgroup_effects = NULL,
      unique_effects = replicate(k, matrix(1, d, d), simplify = FALSE),
      total_effects = replicate(k, matrix(1, d, d), simplify = FALSE),
      tvp_effects = lapply(seq_len(k), function(i) {
        replicate(length(breaks[[i]]), matrix(1, d, d), simplify = FALSE)
      }),
      common_tvp_effects = NULL
    )
  } else if (subgroup) {
    list(
      common_effects = matrix(1, d, d),
      subgroup_effects = replicate(max(subgroup_membership), matrix(1, d, d), simplify = FALSE),
      unique_effects = replicate(k, matrix(1, d, d), simplify = FALSE),
      total_effects = replicate(k, matrix(1, d, d), simplify = FALSE),
      tvp_effects = NULL,
      common_tvp_effects = NULL
    )
  } else {
    list(
      common_effects = matrix(1, d, d),
      subgroup_effects = NULL,
      unique_effects = replicate(k, matrix(1, d, d), simplify = FALSE),
      total_effects = replicate(k, matrix(1, d, d), simplify = FALSE),
      tvp_effects = NULL,
      common_tvp_effects = NULL
    )
  }
}


#' Estimate raw total effects per subject (and per period if TVP)
#'
#' @param Ak List of design matrices per subject
#' @param bk List of response matrices per subject
#' @param weightest Character. Estimation method: "lasso", "alasso", "ridge", or "ols"
#' @param tvp Logical. Whether to estimate time-varying parameters
#' @param breaks List of period indices per subject
#' @param nfolds Numeric. Number of CV folds
#' @param lambda_choice Character. Which lambda to use
#'
#' @return List with total_effects and period_effects (if TVP)
#' @keywords internal
estimate_raw_effects <- function(Ak, bk, weightest, tvp, breaks, nfolds, lambda_choice) {

  # For alasso, run lasso first then refine with adaptive weights
  base_method <- if (weightest == "alasso") "lasso" else weightest

  # Set alpha for glmnet: 1 = lasso, 0 = ridge
  alpha <- switch(base_method,
    "lasso" = 1,
    "ridge" = 0,
    "ols" = NA
  )

  # Build fold structure for CV (only needed for lasso/ridge)
  if (base_method %in% c("lasso", "ridge")) {
    fold_structure <- build_cv_folds(Ak, nfolds)
  } else {
    fold_structure <- NULL
  }

  # Estimate total effects per subject
  total_effects <- lapply(seq_along(Ak), function(g) {
    if (base_method == "ols") {
      estimate_ols(Ak[[g]], bk[[g]])
    } else {
      estimate_glmnet(Ak[[g]], bk[[g]], alpha, fold_structure[[g]], lambda_choice)
    }
  })

  # Refine with adaptive lasso if weightest = "alasso"
  if (weightest == "alasso") {
    total_effects <- lapply(seq_along(Ak), function(g) {
      estimate_glmnet_alasso(
        Ak[[g]], bk[[g]], total_effects[[g]],
        fold_structure[[g]], lambda_choice
      )
    })
  }

  # Estimate per-period effects if TVP
  period_effects <- NULL
  if (tvp) {
    period_effects <- lapply(seq_along(Ak), function(g) {
      lapply(seq_along(breaks[[g]]), function(p) {
        idx <- breaks[[g]][[p]]
        if (base_method == "ols") {
          estimate_ols(Ak[[g]][idx, , drop = FALSE], bk[[g]][idx, , drop = FALSE])
        } else {
          period_folds <- build_period_folds(idx, nfolds)
          # TODO: Change glmnet_intercept to FALSE for consistency with total_effects
          # estimation. Both settings produce equivalent results empirically, but
          # FALSE is more principled (glmnet's intercept is not the right kind for
          # VAR models). Requires regenerating RDS test fixtures.
          # Currently using TRUE to match legacy behavior.
          init_est <- estimate_glmnet(
            Ak[[g]][idx, , drop = FALSE],
            bk[[g]][idx, , drop = FALSE],
            alpha, period_folds, lambda_choice,
            glmnet_intercept = TRUE
          )
          if (weightest == "alasso") {
            estimate_glmnet_alasso(
              Ak[[g]][idx, , drop = FALSE],
              bk[[g]][idx, , drop = FALSE],
              init_est, period_folds, lambda_choice,
              glmnet_intercept = TRUE
            )
          } else {
            init_est
          }
        }
      })
    })
  }

  list(total_effects = total_effects, period_effects = period_effects)
}


#' Build CV fold structure for blocked time series cross-validation
#'
#' @param Ak List of design matrices per subject
#' @param nfolds Number of folds
#'
#' @return List of fold ID vectors per subject
#' @keywords internal
build_cv_folds <- function(Ak, nfolds) {
  make_folds <- function(x, nfolds) {
    split(x, cut(seq_along(x), nfolds, labels = FALSE))
  }

  lapply(Ak, function(A) {
    n <- nrow(A)
    indices <- seq_len(n)
    folds <- make_folds(indices, nfolds)
    # Convert to fold IDs
    unlist(lapply(seq_along(folds), function(fold_num) {
      rep(fold_num, length(folds[[fold_num]]))
    }))
  })
}


#' Build CV folds for a single period
#'
#' @param idx Indices within the period
#' @param nfolds Number of folds
#'
#' @return Fold ID vector
#' @keywords internal
build_period_folds <- function(idx, nfolds) {
  period_length <- length(idx)
  period_indices <- seq_len(period_length)

  make_folds <- function(x, nfolds) {
    split(x, cut(seq_along(x), nfolds, labels = FALSE))
  }

  folds <- make_folds(period_indices, nfolds)
  unlist(lapply(seq_along(folds), function(fold_num) {
    rep(fold_num, length(folds[[fold_num]]))
  }))
}


#' Estimate coefficients using OLS
#'
#' @param A Design matrix (n x d)
#' @param b Response matrix (n x d)
#'
#' @return Coefficient matrix (d x d)
#' @keywords internal
estimate_ols <- function(A, b) {
  # β = (A'A)^{-1} A'b
  # Each column of b is a separate response
  t(solve(t(A) %*% A) %*% t(A) %*% b)
}


#' Estimate coefficients using glmnet (lasso or ridge)
#'
#' @param A Design matrix (n x d)
#' @param b Response matrix (n x d)
#' @param alpha Elastic net mixing parameter (1 = lasso, 0 = ridge)
#' @param folds Fold ID vector for CV
#' @param lambda_choice Which lambda to use ("lambda.min" or "lambda.1se")
#' @param glmnet_intercept Logical. Whether glmnet should fit an intercept (default FALSE)
#'
#' @return Coefficient matrix (d x d)
#' @keywords internal
estimate_glmnet <- function(A, b, alpha, folds, lambda_choice, glmnet_intercept = FALSE) {
  n_responses <- ncol(b)
  n_predictors <- ncol(A)

  fit <- glmnet::cv.glmnet(
    x = diag(n_responses) %x% A,
    y = as.vector(b),
    family = "gaussian",
    alpha = alpha,
    standardize = FALSE,
    intercept = glmnet_intercept,
    foldid = rep(folds, n_responses)
  )

  coefs <- coef(fit, s = lambda_choice)[-1]
  matrix(coefs, n_responses, n_predictors, byrow = TRUE)
}


#' Estimate coefficients using adaptive LASSO via glmnet penalty.factor
#'
#' Two-stage procedure: uses initial LASSO estimates to construct adaptive
#' penalty weights (1/|init|), then runs cv.glmnet with those weights.
#' This produces sparser pilot estimates than plain LASSO, particularly
#' at larger sample sizes where LASSO tends to include too many variables.
#'
#' @param A Design matrix (n x d)
#' @param b Response matrix (n x d)
#' @param init_coefs Initial coefficient matrix from LASSO (d x d)
#' @param folds Fold ID vector for CV
#' @param lambda_choice Which lambda to use ("lambda.min" or "lambda.1se")
#' @param glmnet_intercept Logical. Whether glmnet should fit an intercept (default FALSE)
#' @param adapower Power for adaptive weighting (default 1)
#'
#' @return Coefficient matrix (d x d)
#' @keywords internal
estimate_glmnet_alasso <- function(A, b, init_coefs, folds, lambda_choice,
                                   glmnet_intercept = FALSE, adapower = 1) {
  n_responses <- ncol(b)
  n_predictors <- ncol(A)

  # Build penalty factors from initial estimates
  # init_coefs is d_response x d_predictor (row = response equation)
  # glmnet stacked system has d_response*d_predictor coefficients
  # ordered: all predictors for response 1, then response 2, etc. (byrow=TRUE in matrix())
  penalty_vec <- as.vector(t(1 / abs(init_coefs)^adapower))

  # Cap infinite penalties at a finite value that is large relative to the
  # finite penalties, but not so extreme that glmnet's internal rescaling
  # (penalty.factor is rescaled to sum = nvars) makes the finite penalties
  # effectively zero. A cap of 10x the max finite penalty keeps the ratio
  # meaningful after rescaling.
  finite_vals <- penalty_vec[is.finite(penalty_vec)]
  if (length(finite_vals) > 0 && max(finite_vals) > 0) {
    cap <- max(finite_vals) * 10
  } else {
    cap <- 1e10
  }
  penalty_vec[is.infinite(penalty_vec)] <- cap

  fit <- glmnet::cv.glmnet(
    x = diag(n_responses) %x% A,
    y = as.vector(b),
    family = "gaussian",
    alpha = 1,
    standardize = FALSE,
    intercept = glmnet_intercept,
    penalty.factor = penalty_vec,
    foldid = rep(folds, n_responses)
  )

  coefs <- coef(fit, s = lambda_choice)[-1]
  matrix(coefs, n_responses, n_predictors, byrow = TRUE)
}


#' Decompose raw effects into common/unique/subgroup/tvp structure
#'
#' @param raw List with total_effects and period_effects
#' @param k Number of subjects
#' @param d Number of variables
#' @param subgroup Logical. Whether subgrouping is enabled
#' @param subgroup_membership Vector of subgroup assignments
#' @param tvp Logical. Whether TVP is enabled
#' @param common_effects Logical. Whether to compute common effects
#' @param common_tvp_effects Logical. Whether to compute common TVP effects
#'
#' @return List of decomposed effect matrices
#' @keywords internal
decompose_effects <- function(raw, k, d, subgroup, subgroup_membership,
                               tvp, common_effects, common_tvp_effects) {

  total <- raw$total_effects
  periods <- raw$period_effects

  # --- Common effects ---
  common <- compute_common_effects(total, periods, k, d, tvp, common_effects)

  # --- Subgroup effects ---
  if (subgroup) {
    subgrp <- compute_subgroup_effects(total, subgroup_membership, common, d)
    unique <- lapply(seq_along(total), function(i) {
      total[[i]] - common - subgrp[[subgroup_membership[i]]]
    })
  } else {
    subgrp <- NULL
    if (k == 1) {
      # k=1: unique is zero (all variation is in common or TVP)
      unique <- list(matrix(0, d, d))
    } else {
      unique <- lapply(total, function(t) {
        if (!is.null(common)) t - common else t
      })
    }
  }

  # --- TVP effects ---
  tvp_eff <- NULL
  common_tvp <- NULL

  if (tvp) {
    tvp_decomp <- decompose_tvp_effects(
      periods, k, d, common, common_effects, common_tvp_effects
    )
    tvp_eff <- tvp_decomp$tvp_effects
    common_tvp <- tvp_decomp$common_tvp_effects
  }

  list(
    common_effects = common,
    subgroup_effects = subgrp,
    unique_effects = unique,
    tvp_effects = tvp_eff,
    common_tvp_effects = common_tvp,
    total_effects = total
  )
}


#' Compute common effects as median across subjects (or periods for k=1 TVP)
#'
#' @keywords internal
compute_common_effects <- function(total, periods, k, d, tvp, include_common) {

  if (!include_common) {
    return(NULL)
  }

  if (k == 1 && tvp) {
    # k=1 TVP: common = median across periods
    num_periods <- length(periods[[1]])
    period_array <- array(unlist(periods[[1]]), c(d, d, num_periods))
    apply(period_array, 1:2, median)
  } else if (k > 1) {
    # k>1: common = median across subjects
    total_array <- array(unlist(total), c(d, d, k))
    apply(total_array, 1:2, median)
  } else {
    # k=1 non-TVP: common = total (only one subject)
    total[[1]]
  }
}


#' Compute subgroup effects as median within subgroup minus common
#'
#' @keywords internal
compute_subgroup_effects <- function(total, subgroup_membership, common, d) {

  n_subgroups <- max(subgroup_membership)
  k <- length(total)

  lapply(seq_len(n_subgroups), function(s) {
    members <- which(subgroup_membership == s)
    if (length(members) == 1) {
      # Single member: subgroup effect = total - common
      total[[members]] - common
    } else {
      # Multiple members: median within subgroup, then subtract common
      subgrp_array <- array(unlist(total[members]), c(d, d, length(members)))
      apply(subgrp_array, 1:2, median) - common
    }
  })
}


#' Decompose TVP effects into common_tvp and unique tvp components
#'
#' @keywords internal
decompose_tvp_effects <- function(periods, k, d, common, include_common, include_common_tvp) {

  num_periods <- length(periods[[1]])

  # --- Common TVP effects (k>1 only) ---
  common_tvp <- NULL
  if (include_common_tvp && k > 1) {
    common_tvp <- lapply(seq_len(num_periods), function(p) {
      # Collect period p estimates from all subjects
      period_estimates <- lapply(seq_len(k), function(i) periods[[i]][[p]])
      period_array <- array(unlist(period_estimates), c(d, d, k))
      apply(period_array, 1:2, median)
    })
  }

  # --- TVP effects (deviations) ---
  if (k == 1) {
    # k=1: tvp = period - common
    if (include_common) {
      tvp_eff <- list(
        lapply(periods[[1]], function(p_mat) p_mat - common)
      )
    } else {
      # No common: tvp = period estimates directly
      tvp_eff <- list(periods[[1]])
    }
  } else {
    # k>1: tvp = period - common_tvp (or period - common if no common_tvp)
    if (!is.null(common_tvp)) {
      tvp_eff <- lapply(seq_len(k), function(i) {
        lapply(seq_len(num_periods), function(p) {
          periods[[i]][[p]] - common_tvp[[p]]
        })
      })
    } else {
      # No common TVP: tvp = period - subject's base (unique + common)
      # This ensures total_t = common + unique + tvp_t
      tvp_eff <- lapply(seq_len(k), function(i) {
        # Subject i's base estimate (average across periods)
        subject_periods <- periods[[i]]
        subject_array <- array(unlist(subject_periods), c(d, d, num_periods))
        subject_base <- apply(subject_array, 1:2, median)

        lapply(seq_len(num_periods), function(p) {
          periods[[i]][[p]] - subject_base
        })
      })
    }
  }

  list(tvp_effects = tvp_eff, common_tvp_effects = common_tvp)
}

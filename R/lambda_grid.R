#' Build the lambda1 path (one column per penalty scenario)
#'
#' Construct the sequence (path) of lambda1 values that will be used for cv.
#'
#' In this package, lambda1 is the main sparsity / regularization strength
#' applied to the common effects. Laer, we scaled these into the unique effects
#' via ratios. Importantly, we don't just pick a single lambda1 — we build
#' a path of candidate lambda1 values, from very large (forces almost all penalized
#' coefficients to zero) down to much smaller (less shrinkage).
#'
#' This path is defined by a decreasing set of values, generated from lambda_max 
#' down to lambda_max / depth, so`depth` is used to control that span and the 
#' default value is 1,000.
#'
#' Note: the models have multiple penalty scenarios as alluded to earlier. We
#' call these scenarios, where each scenario encodes a different way to penalize 
#' the common and unique effects. Those scenarios live along the 3rd dimension of 
#' the array W are typically of length length(ratios)*legnth(lambda1), such that
#' 
#'   dim(W) = [ n_outcomes , n_predictors , n_scenarios ]
#'
#' We build one lambda1 path per scenario. That’s why the return value is a
#' matrix with:
#' 
#'   rows   = number of lambda1 values along the path
#'   columns = number of scenarios
#'
#' - cv.multivar() stores this matrix in object@lambda1.
#' - During CV, for each scenario column i and each row r in that column,
#'   the solver fits the model at that lambda1 and scores prediction error.
#'
#' - We also have object@ratios, where each ratio says how strong to
#'   penalize deviation blocks (individual/subgroup/time-varying) relative to
#'   the common effects. Conceptually this means lambda2 = ratio * lambda1.
#'
#' - We use the ingredients:
#'       * lambda1 grid (here)
#'       * ratios       (constructModel)
#'   to construct the lambda2 grid..
#'
#' Some notes on the arguments:
#' 
#' depth      numeric.
#'            Default is 1000.
#'            Controls how wide the lambda1 path is. The top of the path
#'            is lambda_max (the largest penalty we think we need to zero
#'            everything out). The bottom is lambda_max / depth.
#'            Typical depth ~ 1000.
#'
#' nlam       integer.
#'            Default is 30.
#'            Number of lambda1 values along the path (rows of the output).
#'
#' Y          numeric outcome matrix.
#'            Matrix of outcomes.
#'
#' Z          numeric matrix.
#'            Matrix of predictors.
#'
#' W          numeric array of penalty weights.
#'            Shape: [ n_outcomes , n_predictors , n_scenarios ].
#'
#'            Interpretation:
#'            - W[,,i] is the penalty weight layout for scenario i.
#'              Different scenarios correspond to different relative
#'              penalization of global vs subgroup vs individual vs
#'              time-varying effects.
#'
#' tol        numeric scalar.
#'            Default 1e-4.
#'            Tolerance used in the adaptive branch (lamadapt = TRUE) when
#'            binary-searching for the scenario-specific max lambda.
#'            If lamadapt = FALSE, tol is not used.
#'
#' intercept  logical, default FALSE.
#'            Whether the model includes an intercept column. Intercepts are
#'            usually unpenalized. In lamadapt = TRUE mode, we drop the first
#'            column of the fitted coefficient array when checking if
#'            everything is (approximately) zero, so we don't falsely call
#'            the intercept "nonzero".
#'
#' lamadapt   logical, default FALSE.
#'            If FALSE:
#'                We compute one proxy lambda_max from data (max(Y %*% t(Z))),
#'                and then generate the same log-spaced path for every scenario.
#'
#'            If TRUE:
#'                For each scenario i, we refine its own lambda_max using
#'                weighted penalized fits with W[,,i], using a binary search.
#'                That is closer to "adaptive lasso" logic: strong signals get
#'                downweighted and don't need as huge a lambda to kill them;
#'                weak ones get upweighted and do.
#'
#' k          integer or NULL.
#'            Number of individuals. Only needed if lamadapt = TRUE, because
#'            we pass it to the solver (wlasso) inside the scenario-specific
#'            binary search. Ignored when lamadapt = FALSE.
#'
#' Notes on the returned object:
#' 
#' lambda1    matrix, dimension [nlam x n_scenarios].
#'
#'            Column i is the lambda1 path for scenario i.
#'            Row r is a particular lambda1 magnitude along that path
#'
#'            Later:
#'            - cv.multivar() stores this in object@lambda1,
#'              and cross-validates across [row r, column i].
#'            - Combined with object@ratios, you can reconstruct lambda2
#'              via lambda2 = ratio * lambda1 and plot CV error surfaces.
#'
#' @export
lambda_grid <- function(depth = 1000,
                        nlam= 30,
                        Y,
                        Z,
                        W = NULL,
                        tol = 1e-4,
                        intercept = FALSE,
                        lamadapt = FALSE,
                        k = NULL) {


  if (is.null(W)){
    
    n_scenarios <- nlam*nlam
    W <- array(1, dim = c(ncol(Y), ncol(Z), n_scenarios))
    
  } else {
    # W must be 3D: [outcome, predictor, scenario]
    if (!is.array(W) || length(dim(W)) != 3) {
      stop("multivar ERROR: lambda_grid: 
           W must be a 3D array [n_outcomes, n_predictors, n_scenarios].")
    }
    # number of penalty scenarios encoded in W
    n_scenarios <- dim(W)[3]
  }
  
  if (missing(depth) || missing(nlam) || missing(Y) || missing(Z) || missing(W)) {
    stop("multivar ERROR: lambda_grid: depth, nlam, Y, Z, and W are required.")
  }

  # if we want adaptive lambda_max search, we need k (for wlasso call)
  if (lamadapt && is.null(k)) {
    stop("multivar ERROR: lambda_grid: 'k' must be supplied when lamadapt = TRUE.")
  }

  # preallocate output: rows = lambda path length, cols = scenarios
  lambda1 <- matrix(0, nrow = nlam, ncol = n_scenarios)

  # We need a large lambda to start the path. In lasso logic,
  # lambda_max is the smallest lambda that keeps all penalized coefs at 0.
  #
  # Here we approximate that using the max absolute correlation-like term
  # max(Y %*% t(Z)). This is a heuristic for "what scale of lambda basically
  # nukes everything".
  #
  # NOTE: in lamadapt=TRUE mode we will refine this per-scenario.
  lambdah0 <- max(Y %*% t(Z))

  ## ------------------------- main loop over scenarios -----------------
  for (i in seq_len(n_scenarios)) {

    if (!lamadapt) {

      # ---------------------------------------------------------------
      # NON-ADAPTIVE BRANCH
      # We don't try to tailor lambda_max separately for each scenario.
      # We just reuse lambdah0 for all scenarios.
      # ---------------------------------------------------------------
      lambda_hi <- lambdah0

    } else {
      # ---------------------------------------------------------------
      # ADAPTIVE BRANCH
      #
      # Goal:
      #   Find a scenario-specific lambda_max for scenario i, taking into
      #   account the penalty weights in W[,,i].
      #
      # We approximate this with a binary search:
      #   - Start with an interval [lam_lo, lam_hi] = [0, lambdah0]
      #   - Try lambda_try in the middle
      #   - Fit a very penalized weighted-lasso (wlasso) with W[,,i]
      #   - Check: are all penalized coefficients ~0 (within tol)?
      #       yes -> we can lower lam_hi
      #       no  -> we need to raise lam_lo
      #   - Repeat until lam_hi - lam_lo is tiny
      #
      # After convergence, lam_hi is our scenario-specific lambda_max.
      # ---------------------------------------------------------------

      # size info: number of outcomes = dim(W[,,1])[1]
      d_out <- dim(W[,,1])[1]

      # zero coefficient template that matches what wlasso() expects:
      B0 <- array(
        0,
        dim = c(
          d_out,
          dim(W[,,1])[2],
          1
        )
      )

      lam_lo <- 0
      lam_hi <- lambdah0

      # binary search for lambda_max for scenario i
      while ((lam_hi - lam_lo) > 10 * tol) {

        lambda_try <- 0.5 * (lam_hi + lam_lo)

        B_fit <- wlasso(
          B0,
          Z,
          Y,
          W[,,i, drop = FALSE],
          k       = k,
          d       = d_out,
          lambda  = lambda_try,
          eps     = tol,
          intercept = intercept
        )

        # Drop intercept column when checking if all coefficients are zero
        # (intercept is unpenalized, so it shouldn't affect lambda_max search)
        start_col <- if(intercept) 2 else 1
        penalized_block <- B_fit[, start_col:dim(B_fit)[2], , drop = FALSE]

        if (max(abs(penalized_block)) < tol) {
          # still all basically zero at this lambda_try:
          # we can shrink the high end
          lam_hi <- lambda_try
        } else {
          # some coefficient "woke up":
          # we need more penalty to kill it, move lower bound up
          lam_lo <- lambda_try
        }
      }

      lambda_hi <- lam_hi
    }

    # ---------------------------------------------------------------
    # Having chosen a top-of-path lambda for this scenario (`lambda_hi`),
    # build a log-spaced decreasing sequence of length nlam from
    # lambda_hi down to lambda_hi / depth.
    #
    # Row 1  = lambda_hi      (very strong shrinkage)
    # Row n  = lambda_hi/depth (much weaker shrinkage)
    #
    # This mirrors glmnet's lambda path logic.
    # ---------------------------------------------------------------
    lambda1[, i] <- exp(seq(
      from       = log(lambda_hi),
      to         = log(lambda_hi / depth),
      length.out = nlam
    ))
  }

  ## ------------------------- return -------------------------
  # lambda1 is now an nlam x n_scenarios matrix.
  #
  # - Column i is the lambda1 path for scenario i (i.e. for W[,,i]).
  # - Row r is the r-th lambda setting in that path.
  #
  # Downstream:
  #   object@lambda1 <- lambda1
  #   cv.multivar() will try each row/col combo, and in
  #   combination with object@ratios (which define lambda2 = ratio * lambda1),
  #   will select the best-performing penalty configuration.
  return(lambda1)
}

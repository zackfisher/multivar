#' Construct an object of class multivar
#' 
#' @param data List. A list (length = k) of T by d multivariate time series
#' @param lag Numeric. The VAR order. Default is 1.
#' @param horizon Numeric. Desired forecast horizon. Default is 1. ZF Note: Should probably be zero.
#' @param t1 Numeric. Index of time series in which to start cross validation. If NULL, default is floor(nrow(n)/3) where nk is the time series length for individual k.
#' @param t2 Numeric. Index of times series in which to end cross validation. If NULL, default is floor(2*nrow(n)/3) where nk is the time series length for individual k.
#' @param lambda1 Matrix. Regularization parameter grid. Default is NULL (auto-generated).
#' @param nlambda1 Numeric. Number of lambda1 values to search over. Default is 30.
#' @param n_ratios_subgroup Numeric. Number of ratios_subgroup values to search over. Default is 30.
#' @param depth Numeric. Depth of lambda1 grid construction (lambda_min = lambda_max / depth). Default is 10000.
#' @param tol Numeric. Optimization tolerance (default 1e-4).
#' @param window Numeric. Size of rolling window.   
#' @param standardize Logical. Default is true. Whether to standardize the individual data. Note, if intercept = TRUE and standardize = TRUE, the data is scaled but not de-meaned.
#' @param weightest Character. How to estimate initial coefficients for adaptive weights. Default is "lasso". Other options include "alasso" (adaptive LASSO: two-stage procedure using LASSO pilot then adaptive reweighting for sparser estimates), "ridge", "ols", and "multivar". The "multivar" option fits a standard lasso multivar model first to get structured initial estimates. Note: for k=1 TVP models, "multivar" works well with common_effects=TRUE but not with common_effects=FALSE (use "lasso" instead). Only used when lassotype = "adaptive" (ignored for standard LASSO).
#' @param canonical Logical. Default is false. If true, individual datasets are fit to a VAR(1) model.
#' @param threshold Logical. Default is false. If true, and canonical is true, individual transition matrices are thresholded based on significance.
#' @param lassotype Character. Default is "adaptive". Choices are "standard" or "adaptive" lasso.
#' @param intercept Logical. Default is FALSE.
#' @param W Matrix. Default is NULL. 
#' @param ratios_unique Numeric vector. Penalty ratio for unique effects. Default is NULL.
#' @param ratios_subgroup Numeric vector. Penalty ratio for subgroup effects. Default is NULL. 
#' @param ratios_unique_tvp Numeric vector. Default is NULL. 
#' @param cv Character. Default is "rolling" for rolling window cross-validation. "blocked" is also available for blocked folds cross-validation. If "blocked" is selected the nfolds argument should bbe specified.
#' @param nfolds Numeric. The number of folds for use with "blocked" cross-validation.
#' @param lamadapt Logical. Should the lambdas be calculated adaptively. Default is FALSE.
#' @param subgroup_membership Numeric. Vector of subgroup assignments.
#' @param subgroup Logical. Internal argument whether to run subgrouping algorithm.
#' @param B Matrix. Default is NULL.
#' @param pendiag Logical. Logical indicating whether autoregressive parameters should be penalized. Default is TRUE.
#' @param tvp Logical. Default is FALSE.
#' @param inittvpcoefs List.
#' @param breaks List. A list of length K indicating structural breaks in the time series.
#' @param lambda_choice Character. Which lambda to use for initial coefficient estimation: "lambda.min" (default) or "lambda.1se". lambda.min provides better coefficient recovery, especially for small samples and TVP models; lambda.1se may be too conservative, causing all-zero initial estimates.
#' @param common_effects Logical. Whether to include common effects in TVP models. Only applies when tvp = TRUE. Default is TRUE (include common effects). When FALSE, the model becomes Total = Unique + TVP instead of Total = Common + Unique + TVP. This can be useful when you expect no shared dynamics across subjects.
#' @param common_tvp_effects Logical. Whether to include time-varying common effects in TVP models. Default is NULL, which automatically sets to TRUE when tvp = TRUE and FALSE when tvp = FALSE. Only meaningful when tvp = TRUE.
#' @param save_beta Logical. Whether to retain the full beta coefficient array in the cv.multivar result. Default is TRUE. Set to FALSE to reduce memory usage when only the best-model coefficients (in mats) are needed.
#' @param eps Numeric. FISTA convergence tolerance. Default is 1e-3. Smaller values yield more precise solutions but increase computation time.
#' @param warmstart Logical. Whether to use the previous lambda's solution as the starting point for the next lambda in the FISTA solver. Default is TRUE. Reduces computation time with negligible effect on accuracy.
#' @param stopping_crit Character. FISTA convergence criterion. One of "absolute" (default), "relative", or "objective". "absolute" checks max|B_new - B_old| < eps; "relative" normalizes by max|B_old|; "objective" checks relative change in the objective function.
#' @param selection Character. Model selection criterion. \code{"cv"} (default) uses cross-validated MSFE. \code{"ebic"} skips CV folds entirely, fits once on the full data, and selects by Extended BIC. EBIC is much faster and can improve structure recovery when n >> p.
#' @param ebic_gamma Numeric. EBIC tuning parameter, used when \code{selection = "ebic"}. \code{0} gives standard BIC; \code{0.5} (default) is moderate EBIC; \code{1} is most conservative.
#' @param weight_type Character. Adaptive weight function type. \code{"standard"} (default) uses \code{1/|coef|^gamma} with Inf replaced by 1e10. \code{"bounded"} uses \code{1/(1+|coef|/tau)^gamma} where tau is the median of nonzero |coefs|, producing weights in [0,1] with no infinities.
#' @param ncores Numeric. Number of cores for parallel cross-validation. Default is 1.
#' @param max_grid_size Numeric. Maximum number of hyperparameter combinations. If the full grid exceeds this, dimensions are coarsened proportionally. Default is NULL (no limit).
#' @examples
#' 
#' sim  <- multivar_sim(
#'   k = 2,  # individuals
#'   d = 3,  # number of variables
#'   n = 20, # number of timepoints
#'   prop_fill_com = 0.1, # proportion of paths common
#'   prop_fill_ind = 0.1, # proportion of paths unique
#'   lb = 0.1,  # lower bound on coefficient magnitude
#'   ub = 0.9,  # upper bound on coefficient magnitude
#'   sigma = diag(3) # noise
#' )
#' 
#' plot_sim(sim, plot_type = "common")
#' 
#' model <- constructModel(data = sim$data, weightest = "ols")
#'
#' @import Matrix
#' @export
constructModel <- function( data = NULL,
                            lag = 1,
                            horizon = 0,
                            t1 = NULL,
                            t2 = NULL,
                            lambda1 = NULL,
                            nlambda1 = 30,
                            n_ratios_subgroup = 30,
                            depth = NULL,
                            tol = 1e-4,
                            window = 1,
                            standardize = TRUE,
                            weightest = "maity",
                            canonical = FALSE,
                            threshold = FALSE,
                            lassotype = "adaptive",
                            intercept = FALSE,
                            W = NULL,
                            ratios_unique = NULL,
                            ratios_subgroup = NULL,
                            ratios_unique_tvp = NULL,
                            cv = "blocked",
                            nfolds = 10,
                            lamadapt = FALSE,
                            subgroup_membership = NULL,
                            subgroup = FALSE,
                            B = NULL,
                            pendiag = TRUE,
                            tvp = FALSE,
                            inittvpcoefs = list(),
                            breaks = list(),
                            lambda_choice = "lambda.min",
                            common_effects = TRUE,
                            common_tvp_effects = NULL,
                            save_beta = TRUE,
                            ncores = 1,
                            max_grid_size = NULL,
                            eps = 1e-3,
                            warmstart = TRUE,
                            stopping_crit = "absolute",
                            selection = "cv",
                            ebic_gamma = 0.5,
                            weight_type = "standard",
                            maity_opts = list() ){

  #------------------------------------------------------------------
  # basic checks (unchanged)
  #------------------------------------------------------------------
  if (lag != 1){
    stop("multivar ERROR: Currently only lag of order 1 is supported.")
  }

  if (tol < 0 | tol > 1e-1){
    stop("Tolerance must be positive")
  }

  if (!lambda_choice %in% c("lambda.min", "lambda.1se")){
    stop("multivar ERROR: lambda_choice must be either 'lambda.min' or 'lambda.1se'")
  }

  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0)
    stop("multivar ERROR: eps must be a positive numeric scalar.")
  if (!is.logical(warmstart) || length(warmstart) != 1)
    stop("multivar ERROR: warmstart must be TRUE or FALSE.")
  if (!stopping_crit %in% c("absolute", "relative", "objective"))
    stop("multivar ERROR: stopping_crit must be 'absolute', 'relative', or 'objective'.")
  if (!selection %in% c("cv", "ebic"))
    stop("multivar ERROR: selection must be 'cv' or 'ebic'.")
  if (!is.numeric(ebic_gamma) || length(ebic_gamma) != 1 || ebic_gamma < 0)
    stop("multivar ERROR: ebic_gamma must be a non-negative numeric scalar.")
  if (!weight_type %in% c("standard", "bounded"))
    stop("multivar ERROR: weight_type must be 'standard' or 'bounded'.")

  # Set depth default based on model type
  # TVP models need larger depth because adaptive weights can be extreme (up to 1e10)
  # For TVP, nlambda1 determines depth to maintain consistent lambda1 step ratio (~1.6x)
  # Note: depth only affects the lambda1 grid; ratios_unique/ratios_unique_tvp grids are independent
  if (is.null(depth)){
    depth <- 10000
  }

  # Validate common_effects
  if (!is.logical(common_effects) || length(common_effects) != 1){
    stop("multivar ERROR: common_effects must be TRUE or FALSE.")
  }

  # common_effects only works with TVP
  if (!common_effects && !tvp){
    stop("multivar ERROR: common_effects = FALSE requires tvp = TRUE.")
  }

  # common_effects = FALSE not yet supported with subgroups
  if (!common_effects && subgroup){
    stop("multivar ERROR: common_effects = FALSE is not yet supported with subgroups.")
  }

  # Set common_tvp_effects intelligently if not specified
  if (is.null(common_tvp_effects)) {
    common_tvp_effects <- tvp  # TRUE if tvp=TRUE, FALSE otherwise
  }

  # Validate common_tvp_effects
  if (!is.logical(common_tvp_effects) || length(common_tvp_effects) != 1){
    stop("multivar ERROR: common_tvp_effects must be TRUE or FALSE.")
  }

  # common_tvp_effects = FALSE not yet supported with subgroups
  if (!common_tvp_effects && subgroup && tvp){
    stop("multivar ERROR: common_tvp_effects is not yet supported with subgroups.")
  }

  # weightest = "multivar" with TVP requires explicit breaks
  if (weightest == "multivar" && tvp && length(breaks) == 0){
    stop(paste0(
      "multivar ERROR: weightest = 'multivar' with tvp = TRUE requires explicit breaks. ",
      "Please provide structural break points via the 'breaks' argument, ",
      "or use weightest = 'lasso' (the default) instead."
    ))
  }

  # weightest = "maity" restrictions (Maity et al. 2022 MrLasso approach)
  if (weightest == "maity") {
    if (tvp) stop("multivar ERROR: weightest = 'maity' is not yet supported for TVP models.")
    if (subgroup) stop("multivar ERROR: weightest = 'maity' is not yet supported with subgroups.")
    if (length(data) == 1) stop("multivar ERROR: weightest = 'maity' requires k > 1 (multiple subjects).")
  }

  initcoefs <- list()
  
  #------------------------------------------------------------------
  # prep data
  # dat is a list over individuals; each element has $A, $b, $H
  #------------------------------------------------------------------
  dat <- setup_data(data, standardize, lag, horizon, intercept, tvp, breaks)

  #------------------------------------------------------------------
  # Extract data means for intercept recovery (if intercept=TRUE)
  # Intercepts are recovered post-hoc using: c = mean(b) - Phi * mean(A)
  #------------------------------------------------------------------
  if (intercept) {
    data_means <- lapply(dat, function(x) {
      list(
        mean_A = x$mean_A,
        mean_b = x$mean_b,
        sd_A = x$sd_A,
        sd_b = x$sd_b,
        mean_A_periods = x$mean_A_periods,
        mean_b_periods = x$mean_b_periods,
        sd_A_periods = x$sd_A_periods,
        sd_b_periods = x$sd_b_periods
      )
    })
  } else {
    data_means <- list()
  }

  #------------------------------------------------------------------
  # per-individual time series lengths, cols, and CV indices
  #------------------------------------------------------------------
  # ntk: number of usable timepoints (rows of b) per individual k
  ntk <- vapply(dat, function(x) { nrow(x$b) }, numeric(1))
  
  # ndk: number of variables (cols of b) per individual k
  ndk <- vapply(dat, function(x) { ncol(x$b) }, numeric(1))
  
  # Default train/test split indices for each individual
  # (old code used single t1/t2 logic assuming equal n;
  #  now we build vectors t1k/t2k that respect each ntk[k])
  default_t1k <- floor(ntk / 3)
  default_t2k <- ntk
  
  # If the user supplied t1/t2:
  # - if length 1, recycle to all k
  # - if length k, use as-is
  # - else we fall back to defaults
  if (!is.null(t1)) {
    if (length(t1) == length(ntk)) {
      t1k <- as.numeric(t1)
    } else if (length(t1) == 1) {
      t1k <- rep(as.numeric(t1), length(ntk))
    } else {
      stop("multivar ERROR: t1 must be length 1 or length k.")
    }
  } else {
    t1k <- default_t1k
  }
  
  if (!is.null(t2)) {
    if (length(t2) == length(ntk)) {
      t2k <- as.numeric(t2)
    } else if (length(t2) == 1) {
      t2k <- rep(as.numeric(t2), length(ntk))
    } else {
      stop("multivar ERROR: t2 must be length 1 or length k.")
    }
  } else {
    t2k <- default_t2k
  }
  
  #------------------------------------------------------------------
  # construct big design matrix A and stacked outcomes b, H
  #------------------------------------------------------------------
  # helper to build sparse column indices
  getj <- function(mat){
    rep(0:(ncol(mat)-1), each = nrow(mat))
  }
  
  k   <- length(dat)                # number of individuals
  p   <- ncol(dat[[1]]$A)           # number of variables per system (assumed common across individuals)

  # Override common_tvp_effects for k=1 (doesn't make sense for single subject)
  if (k == 1 && common_tvp_effects) {
    common_tvp_effects <- FALSE
  }

  # number of rows per individual's A
  ns  <- vapply(dat, function(item){ nrow(item$A) }, numeric(1))
  cns <- cumsum(ns)
  sr  <- c(1, cns[-length(cns)] + 1) - 1  # 0-based starting row offset for each individual in stacked A
  nz  <- tail(cns, 1)                      # total number of rows across all individuals
  
  is <- js <- xs <- NULL  # will hold i,j,x for sparseMatrix()
  
  #-------------------------------------------------
  # subgroup logic (unchanged except we use ntk[i])
  #-------------------------------------------------
  if (!is.null(subgroup_membership)){
    subgroup <- TRUE
  } else if (subgroup == TRUE & is.null(subgroup_membership)){
    subgroup_membership <- get_subgroups(
      data      = data,
      nlambda1  = nlambda1,
      pendiag   = pendiag
    )
  } else {
    subgroup <- FALSE
    subgroup_membership <- rep(1, k)
  }
  
  if (subgroup){
    # replicate cluster label for each row/column position of subject i
    clust <- lapply(seq_along(subgroup_membership), function(i){
      rep(subgroup_membership[i], ntk[i] * p)
    })
  }
  
  #-------------------------------------------------
  # build block structure for A
  #-------------------------------------------------
  for(ii in 1:k){

    # "group" / common block (skip if common_effects = FALSE and tvp = TRUE)
    if (!(tvp && !common_effects)){
      is <- c(is, sr[ii] + rep(1:nrow(dat[[ii]]$A) - 1L, ncol(dat[[ii]]$A)))
      js <- c(js, getj(dat[[ii]]$A))
      xs <- c(xs, dat[[ii]]$A@x)
    }
    
    if (subgroup){
      # subgroup block
      is <- c(is, sr[ii] + rep(1:nrow(dat[[ii]]$A) - 1L, ncol(dat[[ii]]$A)))
      js <- c(js, getj(dat[[ii]]$A) + clust[[ii]] * p)
      xs <- c(xs, dat[[ii]]$A@x)
      
      # individual-specific block
      is <- c(is, sr[ii] + rep(1:nrow(dat[[ii]]$A) - 1L, ncol(dat[[ii]]$A)))
      js <- c(js, getj(dat[[ii]]$A) + p * (ii + length(unique(subgroup_membership))))
      xs <- c(xs, dat[[ii]]$A@x)
      
    } else {
      # individual-specific block (no subgroup layer)
      is <- c(is, sr[ii] + rep(1:nrow(dat[[ii]]$A) - 1L, ncol(dat[[ii]]$A)))
      # Adjust offset: if no common block, shift indices down by p
      if (tvp && !common_effects){
        js <- c(js, getj(dat[[ii]]$A) + p * (ii - 1))
      } else {
        js <- c(js, getj(dat[[ii]]$A) + p * ii)
      }
      xs <- c(xs, dat[[ii]]$A@x)
    }
  }
  
  Ak <- lapply(dat, "[[", "A")
  bk <- lapply(dat, "[[", "b")
  Hk <- lapply(dat, "[[", "H")
  
  if (subgroup){
    dims_a <- c(nz, p * ((k + 1) + length(unique(subgroup_membership))))
  } else {
    # Reduce dimensions by p if skipping common effects
    if (tvp && !common_effects){
      dims_a <- c(nz, p * k)
    } else {
      dims_a <- c(nz, p * (k + 1))
    }
  }
  
  if (k == 1){
    A <- Matrix(Ak[[1]], sparse = TRUE)
  } else {
    A <- sparseMatrix(
      i      = is,
      j      = js,
      x      = xs,
      index1 = FALSE,
      dims   = dims_a
    )
  }
  
  b <- as.matrix(do.call(rbind, bk))
  H <- as.matrix(do.call(rbind, Hk))
  
  #------------------------------------------------------------------
  # tvp extension (time-varying parameters)
  # this part already naturally respects ntk[i]; we just make sure
  # we split using each subject's own length ntk[i]
  #------------------------------------------------------------------
  if (tvp){

    splitAt <- function(x, pos){
      unname(split(x, cumsum(seq_along(x) %in% pos)))
    }

    # Track whether user provided breaks or we're using defaults
    user_provided_breaks <- (length(breaks) > 0)

    # if user didn't provide breaks, default is each timepoint as a "break"
    # NOTE: This default may not be suitable for TVP estimation (creates single-obs periods)
    # TODO: Consider changing default to create reasonable-sized periods
    if (!user_provided_breaks){
      breaks <- lapply(1:k, function(j){
        seq(from = 1, to = ntk[j], by = 1)
      })
    }

    #------------------------------------------------------------------
    # Validate breaks parameter (only for user-provided breaks)
    #------------------------------------------------------------------
    if (user_provided_breaks){

      # Check breaks is a list of length k
      if (!is.list(breaks)){
        stop("multivar ERROR: breaks must be a list when tvp = TRUE.")
      }

      if (length(breaks) != k){
        stop(sprintf(
          "multivar ERROR: breaks must be a list of length k=%d (got length %d).",
          k, length(breaks)
        ))
      }

      # Validate each subject's breaks before converting to windows
      for (i in seq_along(breaks)){

        # Check breaks[i] is numeric
        if (!is.numeric(breaks[[i]])){
          stop(sprintf(
            "multivar ERROR: breaks[[%d]] must be a numeric vector.",
            i
          ))
        }

        # Check breaks are within valid range
        if (any(breaks[[i]] < 1) || any(breaks[[i]] > ntk[i])){
          stop(sprintf(
            "multivar ERROR: breaks[[%d]] contains indices outside valid range [1, %d].",
            i, ntk[i]
          ))
        }

        # Check breaks are sorted
        if (is.unsorted(breaks[[i]])){
          stop(sprintf(
            "multivar ERROR: breaks[[%d]] must be in chronological order.",
            i
          ))
        }

        # Check for duplicates
        if (any(duplicated(breaks[[i]]))){
          stop(sprintf(
            "multivar ERROR: breaks[[%d]] contains duplicate indices.",
            i
          ))
        }
      }
    }
    
    # convert each subject's break vector into list of contiguous windows
    breaks <- lapply(1:k, function(j){
      splitAt(seq_len(ntk[j]), breaks[[j]])
    })

    #------------------------------------------------------------------
    # Validate period lengths after window conversion (only for user-provided breaks)
    #------------------------------------------------------------------
    if (user_provided_breaks){
      min_period_length <- Inf
      short_period_warnings <- character(0)

      for (i in seq_along(breaks)){
        for (j in seq_along(breaks[[i]])){
          period_len <- length(breaks[[i]][[j]])
          min_period_length <- min(min_period_length, period_len)

          # Check if period is too short for VAR estimation
          # Need at least d+2 observations (d predictors per equation + minimal df)
          min_obs_needed <- ndk[i] + 2
          if (period_len < min_obs_needed){
            stop(sprintf(
              "multivar ERROR: Period %d for subject %d has only %d observations but needs at least %d (d+2 = %d+2) for VAR estimation.",
              j, i, period_len, min_obs_needed, ndk[i]
            ))
          }

          # Warn if period is too short for reliable CV
          if (period_len < nfolds){
            short_period_warnings <- c(short_period_warnings, sprintf(
              "Period %d for subject %d has only %d observations but nfolds=%d. CV will use fewer folds for this period.",
              j, i, period_len, nfolds
            ))
          }
        }
      }

      # Issue all warnings at once
      if (length(short_period_warnings) > 0){
        warning(paste(c(
          "Some periods are shorter than nfolds:",
          short_period_warnings
        ), collapse = "\n  "))
      }

      # Provide guidance for very short periods
      if (min_period_length < 30){
        message(sprintf(
          "Shortest TVP period has %d observations. Consider using fewer CV folds (nfolds < %d) for more stable estimation.",
          min_period_length,
          max(3, floor(min_period_length / 3))
        ))
      }
    }

    # Add time-varying columns to A
    if (k > 1 && common_tvp_effects){
      # k>1 with common TVP: add both common TVP (shared) and unique TVP (block-diagonal)

      # 1. Build common TVP columns (shared across all subjects)
      # Each period gets one set of columns, used by all subjects in that period
      num_periods <- length(breaks[[1]])  # assume same structure for all

      # Create sparse matrix for common TVP
      # We'll build this by stacking each subject's connection to common TVP periods
      common_tvp_parts <- lapply(seq_along(Ak), function(i){
        # For subject i, create columns for each period
        do.call(
          cbind,
          lapply(1:num_periods, function(p){
            # For period p, create columns for each variable
            do.call(
              cbind,
              lapply(1:ncol(Ak[[i]]), function(j){
                # Get rows for subject i in period p
                window <- breaks[[i]][[p]]
                # Create sparse column: 1s in period p rows, 0s elsewhere
                col_data <- numeric(ntk[i])
                col_data[window] <- Ak[[i]][window, j]
                Matrix(col_data, nrow = ntk[i], ncol = 1, sparse = TRUE)
              })
            )
          })
        )
      })

      # Stack all subjects vertically (all use the same common TVP columns)
      common_tvp_block <- do.call(rbind, common_tvp_parts)

      # 2. Build unique TVP columns (block-diagonal, subject-specific)
      unique_tvp_block <- Matrix::bdiag(
        lapply(seq_along(Ak), function(i){
          # for subject i, build block-diagonal "tvp" regressors
          tvp_i <- do.call(
            cbind,
            lapply(1:ncol(Ak[[i]]), function(j){
              Matrix::bdiag(
                lapply(breaks[[i]], function(window){
                  Ak[[i]][window, j, drop = FALSE]
                })
              )
            })
          )
          tvp_i
        })
      )

      # 3. Combine: A + common TVP + unique TVP
      A <- Matrix(
        cbind(A, common_tvp_block, unique_tvp_block),
        sparse = TRUE
      )

    } else {
      # k=1 or common_tvp_effects=FALSE: only unique TVP (block-diagonal by period)
      # Structure: [period1_vars | period2_vars | ...] for each subject
      # Each period block has d columns, with nonzeros only in that period's rows

      # Build TVP block-diagonal columns (period-first, then variable)
      tvp_block <- Matrix::bdiag(
        lapply(seq_along(Ak), function(i){

          # for subject i, build block-diagonal "tvp" regressors
          # Iterate: periods first, then variables within each period
          tvp_i <- do.call(
            cbind,
            lapply(breaks[[i]], function(window){
              # For this period, create d columns (one per variable)
              do.call(cbind, lapply(1:ncol(Ak[[i]]), function(j){
                col_data <- numeric(ntk[i])
                col_data[window] <- Ak[[i]][window, j]
                Matrix(col_data, ncol = 1, sparse = TRUE)
              }))
            })
          )
          tvp_i
        })
      )

      # For k=1 TVP with common_effects=FALSE: skip time-invariant columns
      if (k == 1 && !common_effects) {
        A <- Matrix(tvp_block, sparse = TRUE)
      } else {
        A <- Matrix(cbind(A, tvp_block), sparse = TRUE)
      }
    }
  }
  
  #------------------------------------------------------------------
  # tuning parameter grids (unchanged; depends on k, ntk, subgroup)
  #------------------------------------------------------------------
  if (is.null(lambda1)){

    # ratios_unique: lambda2 / lambda1
    # For k=1 TVP, use number of periods instead of k for ratio computation
    # This makes k=1 TVP behave like k=num_periods non-TVP
    effective_k <- k
    if (tvp && k == 1 && length(breaks) > 0) {
      effective_k <- length(breaks[[1]])
    }

    ratios_unique <- rev(round(
      exp(seq(log(effective_k/sqrt(depth)), log(effective_k*sqrt(depth)), length.out = nlambda1)),
      digits = 10
    ))
    
    # ratios_subgroup: tau / lambda1 (subgroup penalty scaling)
    if (subgroup){
      ratios_subgroup <- rev(round(
        exp(seq(
          log(max(subgroup_membership) / depth),
          log(max(subgroup_membership)),
          length.out = n_ratios_subgroup
        )),
        digits = 10
      ))
    } else {
      ratios_subgroup <- rep(1, n_ratios_subgroup)
    }
    
    # ratios_unique_tvp: scaling for tvp penalty (uses k instead of ntk)
    if (tvp){
      if (k == 1) {
        if (common_effects) {
          # For k=1 TVP with common_effects, ratios_unique_tvp scales common weights
          # Max ratio = P (number of periods), analogous to k>1 where max = k
          # Common is estimated from P periods, so can afford P× less penalization
          n_periods <- length(breaks[[1]])
          nalpha <- nlambda1
          ratios_unique_tvp <- rev(round(
            exp(seq(
              log(n_periods/depth),
              log(n_periods),
              length.out = nalpha
            )),
            digits = 10
          ))
        } else {
          # For k=1 TVP without common_effects, no common weights to scale
          nalpha <- 1
          ratios_unique_tvp <- 1
        }
      } else {
        nalpha <- nlambda1
        ratios_unique_tvp <- rev(round(
          exp(seq(
            log(k/depth),
            log(k),
            length.out = nalpha
          )),
          digits = 10
        ))
      }
    } else {
      nalpha <- nlambda1
      ratios_unique_tvp <- rep(1, nalpha)
    }

    # ratios_common_tvp: scaling for common TVP penalty
    if (tvp && common_tvp_effects){
      nbeta <- nlambda1
      ratios_common_tvp <- rev(round(
        exp(seq(
          log(k/depth),
          log(k),
          length.out = nbeta
        )),
        digits = 10
      ))
    } else {
      nbeta <- nlambda1
      ratios_common_tvp <- rep(1, nbeta)
    }

    # Coarsen grids if total size exceeds max_grid_size
    if (!is.null(max_grid_size) && max_grid_size > 0) {
      grid_dims <- c(
        lambda = nlambda1,
        ratio = length(ratios_unique),
        ratio_tvp = length(ratios_unique_tvp),
        ratio_subgroup = length(ratios_subgroup),
        ratio_common_tvp = length(ratios_common_tvp)
      )
      # Only count dimensions > 1 as "active"
      active <- grid_dims[grid_dims > 1]
      total_grid <- prod(grid_dims)

      if (total_grid > max_grid_size) {
        # Proportionally reduce each active dimension
        reduction_factor <- (max_grid_size / total_grid)^(1 / length(active))
        thin_grid <- function(x, n_new) {
          if (length(x) <= n_new) return(x)
          idx <- round(seq(1, length(x), length.out = n_new))
          x[idx]
        }

        if (nlambda1 > 1) {
          nlambda1 <- max(2L, as.integer(ceiling(nlambda1 * reduction_factor)))
        }
        if (length(ratios_unique) > 1) {
          ratios_unique <- thin_grid(ratios_unique,
            max(2L, as.integer(ceiling(length(ratios_unique) * reduction_factor))))
        }
        if (length(ratios_unique_tvp) > 1) {
          ratios_unique_tvp <- thin_grid(ratios_unique_tvp,
            max(2L, as.integer(ceiling(length(ratios_unique_tvp) * reduction_factor))))
        }
        if (length(ratios_subgroup) > 1) {
          ratios_subgroup <- thin_grid(ratios_subgroup,
            max(2L, as.integer(ceiling(length(ratios_subgroup) * reduction_factor))))
        }
        if (length(ratios_common_tvp) > 1) {
          ratios_common_tvp <- thin_grid(ratios_common_tvp,
            max(2L, as.integer(ceiling(length(ratios_common_tvp) * reduction_factor))))
        }

        message(sprintf("Grid coarsened from %d to %d combinations (max_grid_size=%d)",
                        total_grid,
                        nlambda1 * length(ratios_unique) * length(ratios_unique_tvp) *
                          length(ratios_subgroup) * length(ratios_common_tvp),
                        max_grid_size))
      }
    }

    lambda1 <- matrix(0, nlambda1, length(ratios_unique))

  } else {
    # User provided explicit lambda1 values
    nlambda1 <- length(lambda1)
    lambda1 <- matrix(lambda1, nrow = 1)

    # Use user-provided ratios if specified, otherwise default to 1
    ratios_unique      <- if (is.null(ratios_unique)) 1 else ratios_unique
    ratios_unique_tvp  <- if (is.null(ratios_unique_tvp)) 1 else ratios_unique_tvp
    ratios_subgroup    <- if (is.null(ratios_subgroup)) 1 else ratios_subgroup
    ratios_common_tvp  <- 1
  }
  
  #------------------------------------------------------------------
  # W: weight matrix for penalties (unchanged)
  #------------------------------------------------------------------
  W <- matrix(1, nrow = ncol(bk[[1]]), ncol = ncol(A))
  
  #------------------------------------------------------------------
  # B: container for fitted coefficients across tuning params
  #------------------------------------------------------------------
  if (!subgroup){
    if (k == 1){
      B <- array(
        0,
        dim = c(
          ndk[1],
          p,
          nlambda1 * length(ratios_unique)
        )
      )
    } else {
      B <- array(
        0,
        dim = c(
          ndk[1],
          p * (k + 1),
          nlambda1 * length(ratios_unique)
        )
      )
    }
  } else {
    B <- array(
      0,
      dim = c(
        ndk[1],
        p * (k + max(subgroup_membership) + 1),
        nlambda1 * length(ratios_unique) * length(ratios_subgroup)
      )
    )
  }
  
  
  if (tvp){

    # this is not tested intercept coverage


    # overwrite B to include tvp columns.
    # additional predictors for subject i:
    #   ncol(Ak[[i]]) * length(breaks[[i]])
    # total extra cols = sum_i ncol(Ak[[i]]) * length(breaks[[i]])
    extra_cols <- sum(sapply(seq_along(Ak), function(i){
      ncol(Ak[[i]]) * length(breaks[[i]])
    }))

    # For k=1, no separate common/unique effects needed
    # For k>1, use standard (k+1) structure for common + unique
    if (k == 1){
      # k=1 TVP with common_effects=FALSE: no time-invariant columns
      if (common_effects){
        base_cols <- ndk[1]  # time-invariant columns
      } else {
        base_cols <- 0  # no time-invariant columns, only TVP
      }
      # k=1: no ratios_unique needed (no common/unique distinction)
      # For TVP k=1, no common_tvp_effects (no multi-subject decomposition)
      grid_size <- nlambda1 * length(ratios_unique_tvp)
      common_tvp_cols <- 0  # k=1 has no common/unique TVP distinction
      unique_tvp_cols <- extra_cols
    } else {
      # k>1: need to account for common_tvp_effects
      # Base columns (common + unique base effects)
      if (common_effects){
        base_cols <- ndk[1] * (k + 1)
      } else {
        base_cols <- ndk[1] * k
      }

      # TVP columns
      if (common_tvp_effects){
        # Common TVP: one set of columns per period (shared across subjects)
        num_periods <- length(breaks[[1]])  # assume same structure for all subjects
        common_tvp_cols <- num_periods * ndk[1]
        # Unique TVP: per-subject, per-period (block-diagonal)
        unique_tvp_cols <- extra_cols
      } else {
        # No common TVP: all TVP is unique (current behavior)
        common_tvp_cols <- 0
        unique_tvp_cols <- extra_cols
      }

      # Grid size calculation
      if (common_effects && common_tvp_effects){
        # 4-layer: common base, unique base, common TVP, unique TVP
        grid_size <- nlambda1 * length(ratios_unique) * length(ratios_unique_tvp) * length(ratios_common_tvp)
      } else if (common_effects && !common_tvp_effects){
        # 3-layer: common base, unique base, unique TVP (current default)
        grid_size <- nlambda1 * length(ratios_unique) * length(ratios_unique_tvp)
      } else if (!common_effects && common_tvp_effects){
        # 3-layer: unique base, common TVP, unique TVP
        grid_size <- nlambda1 * length(ratios_unique_tvp) * length(ratios_common_tvp)
      } else {
        # 2-layer: unique base, unique TVP (current common_effects=FALSE)
        grid_size <- nlambda1 * length(ratios_unique_tvp)
      }
    }

    B <- array(
      0,
      dim = c(
        ndk[1],
        base_cols + common_tvp_cols + unique_tvp_cols,
        grid_size
      )
    )
  }
  
  #------------------------------------------------------------------
  # Build matrix specification (single source of truth for column indices)
  #------------------------------------------------------------------
  spec <- build_matrix_spec(
    k = k,
    d = ndk[1],
    n = ntk,
    tvp = tvp,
    common_effects = common_effects,
    subgroup = subgroup,
    common_tvp_effects = common_tvp_effects,
    breaks = if (tvp) breaks else NULL,
    subgroup_membership = if (subgroup) subgroup_membership else NULL
  )

  # Verify spec matches actual A matrix dimensions
  if (spec$cols$total != ncol(A)) {
    stop(sprintf(
      "multivar ERROR: matrix spec mismatch. spec$cols$total=%d but ncol(A)=%d",
      spec$cols$total, ncol(A)
    ))
  }
  if (spec$rows$total != nrow(A)) {
    stop(sprintf(
      "multivar ERROR: matrix spec mismatch. spec$rows$total=%d but nrow(A)=%d",
      spec$rows$total, nrow(A)
    ))
  }

  #------------------------------------------------------------------
  # build and return S4 multivar object
  #------------------------------------------------------------------
  obj <- new("multivar",
             k  = k,
             n  = ntk,
             d  = ndk[1],
             Ak = Ak,
             bk = bk,
             Hk = Hk,
             A  = as.matrix(A),
             b  = b,
             H  = H,
             horizon = horizon,
             t1 = t1k,
             t2 = t2k,
             t1k = t1k,
             t2k = t2k,
             ntk = ntk,
             ndk = ndk,
             lambda1 = lambda1,
             nlambda1 = nlambda1,
             n_ratios_subgroup = n_ratios_subgroup,
             depth = depth,
             tol = tol,
             window = window,
             weightest = weightest,
             canonical  = canonical,
             threshold  = threshold,
             lassotype = lassotype,
             intercept = intercept,
             data_means = data_means,
             W = W,
             ratios_unique = ratios_unique,
             ratios_subgroup = ratios_subgroup,
             ratios_unique_tvp = ratios_unique_tvp,
             ratios_common_tvp = ratios_common_tvp,
             cv = cv,
             nfolds = nfolds,
             lamadapt = lamadapt,
             subgroup_membership = subgroup_membership,
             subgroup = subgroup,
             B = B,
             initcoefs = initcoefs,
             pendiag = pendiag,
             tvp = tvp,
             inittvpcoefs = inittvpcoefs,
             breaks = breaks,
             lambda_choice = lambda_choice,
             common_effects = common_effects,
             common_tvp_effects = common_tvp_effects,
             save_beta = save_beta,
             ncores = ncores,
             spec = unclass(spec),
             eps = eps,
             warmstart = warmstart,
             stopping_crit = stopping_crit,
             selection = selection,
             ebic_gamma = ebic_gamma,
             weight_type = weight_type,
             maity_opts = maity_opts
  )

  return(obj)
}

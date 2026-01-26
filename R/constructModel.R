#' Construct an object of class multivar
#' 
#' @param data List. A list (length = k) of T by d multivariate time series
#' @param lag Numeric. The VAR order. Default is 1.
#' @param horizon Numeric. Desired forecast horizon. Default is 1. ZF Note: Should probably be zero.
#' @param t1 Numeric. Index of time series in which to start cross validation. If NULL, default is floor(nrow(n)/3) where nk is the time series length for individual k.
#' @param t2 Numeric. Index of times series in which to end cross validation. If NULL, default is floor(2*nrow(n)/3) where nk is the time series length for individual k.
#' @param lambda1 Matrix. Regularization parameter 1. Default is NULL.
#' @param lambda2 Matrix. Regularization parameter 2. Default is NULL.
#' @param tau Matrix. Regularization parameter for subgroup effects.
#' @param nlambda1 Numeric. Number of lambda1 values to search over. Default is 30.
#' @param nlambda2 Numeric. Number of lambda2 values to search over. Default is 30.
#' @param ntau Numeric. Number of tau values to search over. Default is 30.
#' @param depth Numeric. Depth of grid construction. Default is 1000.
#' @param tol Numeric. Optimization tolerance (default 1e-4).
#' @param window Numeric. Size of rolling window.   
#' @param standardize Logical. Default is true. Whether to standardize the individual data. Note, if intercept = TRUE and standardize = TRUE, the data is scaled but not de-meaned.
#' @param weightest Character. How to estimate initial coefficients for adaptive weights. Default is "lasso". Other options include "ridge" and "ols". Only used when lassotype = "adaptive" (ignored for standard LASSO).
#' @param canonical Logical. Default is false. If true, individual datasets are fit to a VAR(1) model.
#' @param threshold Logical. Default is false. If true, and canonical is true, individual transition matrices are thresholded based on significance.
#' @param lassotype Character. Default is "adaptive". Choices are "standard" or "adaptive" lasso.
#' @param intercept Logical. Default is FALSE.
#' @param pen_common_intercept Logical. Default is FALSE. Whether to penalize the common intercept.
#' @param pen_unique_intercept Logical. Default is TRUE. Whether to penalize individual-specific intercept deviations.
#' @param W Matrix. Default is NULL. 
#' @param ratios Numeric vector. Default is NULL. 
#' @param ratiostau Numeric vector. Default is NULL. 
#' @param ratiosalpha Numeric vector. Default is NULL. 
#' @param cv Character. Default is "rolling" for rolling window cross-validation. "blocked" is also available for blocked folds cross-validation. If "blocked" is selected the nfolds argument should bbe specified.
#' @param nfolds Numeric. The number of folds for use with "blocked" cross-validation.
#' @param thresh Numeric. Post-estimation threshold for setting the individual-level coefficients to zero if their absolute value is smaller than the value provided. Default is zero.
#' @param lamadapt Logical. Should the lambdas be calculated adaptively. Default is FALSE.
#' @param subgroup_membership Numeric. Vector of subgroup assignments.
#' @param subgroup Logical. Internal argument whether to run subgrouping algorithm.
#' @param B Matrix. Default is NULL.
#' @param pendiag Logical. Logical indicating whether autoregressive parameters should be penalized. Default is TRUE.
#' @param tvp Logical. Default is FALSE.
#' @param inittvpcoefs List.
#' @param breaks List. A list of length K indicating structural breaks in the time series.
#' @param lambda_choice Character. Which lambda to use for initial coefficient estimation: "lambda.1se" (default) or "lambda.min". lambda.1se is recommended for sparser, more stable solutions; lambda.min may overfit.
#' @param common_effects Logical. Whether to include common effects in TVP models. Only applies when tvp = TRUE. Default is TRUE (include common effects). When FALSE, the model becomes Total = Unique + TVP instead of Total = Common + Unique + TVP. This can be useful when you expect no shared dynamics across subjects.
#' @param common_tvp_effects Logical. Whether to include time-varying common effects in TVP models. Default is NULL, which automatically sets to TRUE when tvp = TRUE and FALSE when tvp = FALSE. Only meaningful when tvp = TRUE.
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
                            lambda2 = NULL,
                            tau = NULL,
                            nlambda1 = 30,
                            nlambda2 = 30,
                            ntau = 30,
                            depth = 1000,
                            tol = 1e-4,
                            window = 1,
                            standardize = TRUE,
                            weightest = "lasso",
                            canonical = FALSE,
                            threshold = FALSE,
                            lassotype = "adaptive",
                            intercept = FALSE,
                            pen_common_intercept = FALSE,
                            pen_unique_intercept = TRUE,
                            W = NULL,
                            ratios = NULL,
                            ratiostau = NULL,
                            ratiosalpha = NULL,
                            cv = "blocked",
                            nfolds = 10,
                            thresh = 0,
                            lamadapt = FALSE,
                            subgroup_membership = NULL,
                            subgroup = FALSE,
                            B = NULL,
                            pendiag = TRUE,
                            tvp = FALSE,
                            inittvpcoefs = list(),
                            breaks = list(),
                            lambda_choice = "lambda.1se",
                            common_effects = TRUE,
                            common_tvp_effects = NULL ){

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

  initcoefs <- list()
  
  #------------------------------------------------------------------
  # prep data
  # dat is a list over individuals; each element has $A, $b, $H
  #------------------------------------------------------------------
  dat <- setup_data(data, standardize, lag, horizon, intercept)
  
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
      nlambda2  = nlambda2,
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
          # Need at least d*(d+1) observations to estimate d*d transition matrix + d intercepts
          min_obs_needed <- ndk[i] * (ndk[i] + 1)
          if (period_len < min_obs_needed){
            stop(sprintf(
              "multivar ERROR: Period %d for subject %d has only %d observations but needs at least %d (d*(d+1) = %d*%d) for VAR estimation.",
              j, i, period_len, min_obs_needed, ndk[i], ndk[i] + 1
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
      # k=1 or common_tvp_effects=FALSE: only unique TVP (block-diagonal)
      # This is the current behavior
      A <- Matrix(
        cbind(
          A,
          Matrix::bdiag(
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
        ),
        sparse = TRUE
      )
    }
  }
  
  #------------------------------------------------------------------
  # tuning parameter grids (unchanged; depends on k, ntk, subgroup)
  #------------------------------------------------------------------
  if (is.null(lambda1) & is.null(lambda2)){
    
    # ratios: lambda2 / lambda1
    ratios <- rev(round(
      exp(seq(log(k/depth), log(k), length.out = nlambda1)),
      digits = 10
    ))
    
    # ratiostau: tau / lambda1 (subgroup penalty scaling)
    if (subgroup){
      ratiostau <- rev(round(
        exp(seq(
          log(max(subgroup_membership) / depth),
          log(max(subgroup_membership)),
          length.out = ntau
        )),
        digits = 10
      ))
    } else {
      ratiostau <- rep(1, ntau)
    }
    
    # ratiosalpha: scaling for tvp penalty (uses k instead of ntk)
    if (tvp){
      nalpha <- nlambda1
      ratiosalpha <- rev(round(
        exp(seq(
          log(k/depth),
          log(k),
          length.out = nalpha
        )),
        digits = 10
      ))
    } else {
      nalpha <- nlambda1
      ratiosalpha <- rep(1, nalpha)
    }

    # ratiosbeta: scaling for common TVP penalty
    if (tvp && common_tvp_effects){
      nbeta <- nlambda1
      ratiosbeta <- rev(round(
        exp(seq(
          log(k/depth),
          log(k),
          length.out = nbeta
        )),
        digits = 10
      ))
    } else {
      nbeta <- nlambda1
      ratiosbeta <- rep(1, nbeta)
    }

    lambda1 <- matrix(0, nlambda1, length(ratios))
    lambda2 <- matrix(0, nlambda2, length(ratios))
    tau     <- matrix(0, ntau,     length(ratiostau))
    
  } else {
    
    nlambda1 <- length(lambda1)
    nlambda2 <- length(lambda2)
    
    lambda1 <- matrix(lambda1, nrow = 1)
    lambda2 <- matrix(lambda2, nrow = 1)
    
    ratios       <- c(0)
    tau          <- c(0)
    ratiosalpha  <- if (is.null(ratiosalpha)) 1 else ratiosalpha
    ratiostau    <- if (is.null(ratiostau))   1 else ratiostau
  }
  
  #------------------------------------------------------------------
  # W: weight matrix for penalties (unchanged)
  #------------------------------------------------------------------
  W <- matrix(1, nrow = ncol(bk[[1]]), ncol = ncol(A))
  
  #------------------------------------------------------------------
  # B: container for fitted coefficients across tuning params
  # (The logic here is unchanged; depends on k, subgroup, tvp)
  #------------------------------------------------------------------
  # if (!subgroup){
  #   if (k == 1){
  #     B <- array(
  #       0,
  #       dim = c(
  #         ndk[1],
  #         ndk[1]  + 1,
  #         nlambda1 * length(ratios)
  #       )
  #     )
  #   } else {
  #     B <- array(
  #       0,
  #       dim = c(
  #         ndk[1],
  #         ndk[1] * (k + 1) + 1,
  #         nlambda1 * length(ratios)
  #       )
  #     )
  #   }
  # } else {
  #   B <- array(
  #     0,
  #     dim = c(
  #       ndk[1],
  #       ndk[1] * (k + max(subgroup_membership) + 1) + 1,
  #       nlambda1 * length(ratios) * length(ratiostau)
  #     )
  #   )
  # }
  # B dimensions: rows = outcome vars (d), cols = predictor coefficients

  # p accounts for intercept column when intercept=TRUE
  if (!subgroup){
    if (k == 1){
      B <- array(
        0,
        dim = c(
          ndk[1],
          p,
          nlambda1 * length(ratios)
        )
      )
    } else {
      B <- array(
        0,
        dim = c(
          ndk[1],
          p * (k + 1),
          nlambda1 * length(ratios)
        )
      )
    }
  } else {
    B <- array(
      0,
      dim = c(
        ndk[1],
        p * (k + max(subgroup_membership) + 1),
        nlambda1 * length(ratios) * length(ratiostau)
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
      base_cols <- ndk[1]
      # k=1: no ratios needed (no common/unique distinction)
      # For TVP k=1, no common_tvp_effects (no multi-subject decomposition)
      grid_size <- nlambda1 * length(ratiosalpha)
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
        grid_size <- nlambda1 * length(ratios) * length(ratiosalpha) * length(ratiosbeta)
      } else if (common_effects && !common_tvp_effects){
        # 3-layer: common base, unique base, unique TVP (current default)
        grid_size <- nlambda1 * length(ratios) * length(ratiosalpha)
      } else if (!common_effects && common_tvp_effects){
        # 3-layer: unique base, common TVP, unique TVP
        grid_size <- nlambda1 * length(ratiosalpha) * length(ratiosbeta)
      } else {
        # 2-layer: unique base, unique TVP (current common_effects=FALSE)
        grid_size <- nlambda1 * length(ratiosalpha)
      }
    }

    B <- array(
      0,
      dim = c(
        ndk[1],
        base_cols + common_tvp_cols + unique_tvp_cols + if(intercept) 1 else 0,
        grid_size
      )
    )
  }
  
  #------------------------------------------------------------------
  # build and return S4 multivar object
  #------------------------------------------------------------------
  obj <- new("multivar",
             k  = k,
             n  = ntk,
             d  = ndk,
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
             lambda2 = lambda2,
             tau = tau,
             nlambda1 = nlambda1,
             nlambda2 = nlambda2,
             ntau = ntau,
             depth = depth,
             tol = tol,
             window = window,
             weightest = weightest,
             canonical  = canonical,
             threshold  = threshold,
             lassotype = lassotype,
             intercept = intercept,
             pen_common_intercept = pen_common_intercept,
             pen_unique_intercept = pen_unique_intercept,
             W = W,
             ratios = ratios,
             ratiostau = ratiostau,
             ratiosalpha = ratiosalpha,
             ratiosbeta = ratiosbeta,
             cv = cv,
             nfolds = nfolds,
             thresh = thresh,
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
             common_tvp_effects = common_tvp_effects
  )

  return(obj)
}

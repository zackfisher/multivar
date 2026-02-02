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
#' @param weightest Character. How to estimate initial coefficients for adaptive weights: "lasso", "ridge", "ols", or "multivar". The "multivar" option fits a standard lasso multivar model first to get structured initial estimates. Note: for k=1 TVP models, use "lasso" (not "multivar") for best results. Only used when lassotype = "adaptive" (ignored for standard LASSO).
#' @param subgroup_membership Numeric vector. Vector of subgroup assignments for each subject.
#' @param subgroup Logical. Whether to run subgrouping algorithm.
#' @param nlambda1 Numeric. Number of lambda1 values for the main estimation grid. Not used when developing adaptive weights.
#' @param nlambda2 Numeric. Number of lambda2 values for the main estimation grid. Not used when developing adaptive weights.
#' @param tvp Logical. Whether to estimate time-varying parameters.
#' @param breaks List. A list of length k indicating structural breaks in the time series.
#' @param intercept Logical. Whether to include intercepts in the model.
#' @param nfolds Numeric. The number of folds for cross-validation.
#' @param lambda_choice Character. Which lambda to use from cv.glmnet: "lambda.min"
#'   or "lambda.1se". Passed from constructModel; see constructModel documentation
#'   for usage guidance.
#' @param common_effects Logical. Whether to include common effects in TVP models.
#'   Passed from constructModel.
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
    nlambda2,
    tvp,
    breaks,
    intercept,
    nfolds,
    lambda_choice = "lambda.min",
    common_effects = TRUE,
    common_tvp_effects = TRUE){
  
  # Ak<-object@Ak
  # bk<- object@bk
  # ratios<-object@ratios
  # d<-object@d
  # k<-object@k
  # lassotype<-object@lassotype
  # weightest<-object@weightest
  # subgroup <- object@subgroup
  # subgroupflag <- object@subgroup
  # nlambda1 <- object@nlambda1
  # nlambda2 <- object@nlambda2
  # tvp <- object@tvp
  # breaks <- object@breaks
  # intercept <- object@intercept
  # nfolds <- object@nfolds

  if (lassotype == "standard"){
    
    if(!subgroup){
      res <- list(
        common_effects = matrix(1, d, d),
        subgroup_effects = NULL,
        unique_effects = replicate(k,  matrix(1, d, d)),
        total_effects = replicate(k,  matrix(1, d, d)),
        tvp_effects = NULL
      )
    } else {
      res <- list(
        common_effects = matrix(1, d, d),
        subgroup_effects = replicate(max(subgroup_membership),  matrix(1, d, d)),
        unique_effects = replicate(k,  matrix(1, d, d)),
        total_effects = replicate(k,  matrix(1, d, d)),
        tvp_effects = NULL
      )
    }
    
    if(tvp){
      res <- list(
        common_effects = matrix(1, d, d),
        subgroup_effects = NULL,
        unique_effects = replicate(k,  matrix(1, d, d)),
        total_effects = replicate(k,  matrix(1, d, d)),
        tvp_effects = replicate(k,  matrix(1, nrow(Ak[[1]]), d))
      )
    }
    
    
  } else {
    
    #---------------------------------------------#
    # 1. estimate total effects
    #---------------------------------------------#
    if (weightest == "lasso" | weightest == "ridge") {
      
      make_folds  <- function(x,nfolds) split(x, cut(seq_along(x), nfolds, labels = FALSE))
      final_tmpt <- cumsum(unlist(lapply(Ak,function(x){nrow(x)})))
      first_tmpt <- c(1,(final_tmpt[-length(final_tmpt)]+1))
      subj_indx_list <- lapply(seq_along(first_tmpt), function(g){first_tmpt[g]:final_tmpt[g]})
      cv_list <- lapply(subj_indx_list,function(g){make_folds(g,nfolds)})
      
      glmnet_folds <- lapply(seq_along(cv_list), function(jj){
        lapply(seq_along(cv_list[[jj]]), function(kk){
          rep(kk, length(cv_list[[jj]][[kk]]))
          })
      })
      
    }
    
    if (weightest == "ols") {
        
      total_effects <- lapply(seq_along(Ak),function(g){
        if(intercept){
          
          fit<- vars::VAR(as.matrix(bk[[g]]), p=1, type="const")$varresult
          m <- as.matrix(do.call("rbind",lapply(seq_along(colnames(bk[[g]])), function(x) {
            fit[[x]]$coefficients
          })))
          m <- cbind(m[,ncol(m),drop =F],m[,-ncol(m),drop =F])
          m
          
        } else {
          
          fit<- vars::VAR(as.matrix(Ak[[g]]), p=1, type="none")$varresult
          m <- as.matrix(do.call("rbind",lapply(seq_along(colnames(Ak[[g]])), function(x) {
            fit[[x]]$coefficients
          })))
          m
          
        }
      })
        
    } else if (weightest == "lasso") {
      
      total_effects <- lapply(seq_along(Ak),function(g){
        
        lasso.fit <- glmnet::cv.glmnet(
          diag(d[1]) %x% Ak[[g]], 
          as.vector(bk[[g]]), 
          family = "gaussian", 
          alpha = 1, 
          standardize = FALSE, 
          nfolds = nfolds,
          intercept = FALSE,
          foldid = rep(unlist(glmnet_folds[[g]]), d[1])
        )

        # Use specified lambda (default: lambda.min for better parameter recovery)
        g_coefs <- coef(lasso.fit, s = lambda_choice)[-1]
        
        mat <- matrix(g_coefs, ncol(bk[[g]]), ncol(Ak[[g]]), byrow=TRUE)
        
        mat
        
      })
      
      
      # TVP initial coefficient estimation
      # n_responses = number of response equations (d)
      # n_predictors = number of predictors (d without intercept, d+1 with intercept)
      if(tvp){

        inittvpcoefs <- lapply(seq_along(Ak),function(g){

          n_responses <- ncol(bk[[g]])
          n_predictors <- ncol(Ak[[g]])

          lapply(seq_along(breaks[[g]]), function(q){

            # Create blocked folds for this period
            period_length <- length(breaks[[g]][[q]])
            period_indices <- seq_len(period_length)

            # Use blocked folds for time series data
            period_folds <- make_folds(period_indices, nfolds)

            # Convert to fold IDs for glmnet
            period_foldid <- unlist(lapply(seq_along(period_folds), function(fold_num){
              rep(fold_num, length(period_folds[[fold_num]]))
            }))

            # Replicate fold IDs for each outcome variable (Kronecker structure)
            foldid_glmnet <- rep(period_foldid, n_responses)

            lasso.fit <- glmnet::cv.glmnet(
              diag(n_responses) %x% Ak[[g]][breaks[[g]][[q]],],
              as.vector(bk[[g]][breaks[[g]][[q]],]),
              family = "gaussian",
              alpha = 1,
              standardize = FALSE,
              foldid = foldid_glmnet
            )
            matrix(as.vector(predict(
              lasso.fit,
              diag(n_responses) %x% Ak[[g]][breaks[[g]][[q]],],
              type="coefficients",
              s=lambda_choice
            ))[-1], n_responses, n_predictors, byrow=TRUE)
          })
        })

      }
      
    } else if (weightest == "ridge") {

      total_effects <- lapply(seq_along(Ak),function(g){

        ridge.fit <- glmnet::cv.glmnet(
          diag(d[1]) %x% Ak[[g]],
          as.vector(bk[[g]]),
          family = "gaussian",
          alpha = 0,
          standardize = FALSE,
          nfolds = nfolds,
          intercept = FALSE,
          foldid =  rep(unlist(glmnet_folds[[g]]), d[1])
        )

        g_coefs <- coef(ridge.fit, s = lambda_choice)[-1]

        mat <- matrix(g_coefs, ncol(bk[[g]]), ncol(Ak[[g]]), byrow=TRUE)

        mat

      })

      # TVP initial estimates with ridge
      # n_responses = number of response equations (d)
      # n_predictors = number of predictors (d without intercept, d+1 with intercept)
      if(tvp){

        inittvpcoefs <- lapply(seq_along(Ak),function(g){

          n_responses <- ncol(bk[[g]])
          n_predictors <- ncol(Ak[[g]])

          lapply(seq_along(breaks[[g]]), function(q){

            # Create blocked folds for this period
            period_length <- length(breaks[[g]][[q]])
            period_indices <- seq_len(period_length)

            # Use blocked folds for time series data
            period_folds <- make_folds(period_indices, nfolds)

            # Convert to fold IDs for glmnet
            period_foldid <- unlist(lapply(seq_along(period_folds), function(fold_num){
              rep(fold_num, length(period_folds[[fold_num]]))
            }))

            # Replicate fold IDs for each outcome variable (Kronecker structure)
            foldid_glmnet <- rep(period_foldid, n_responses)

            ridge.fit <- glmnet::cv.glmnet(
              diag(n_responses) %x% Ak[[g]][breaks[[g]][[q]],],
              as.vector(bk[[g]][breaks[[g]][[q]],]),
              family = "gaussian",
              alpha = 0,
              standardize = FALSE,
              foldid = foldid_glmnet
            )
            matrix(as.vector(predict(
              ridge.fit,
              diag(n_responses) %x% Ak[[g]][breaks[[g]][[q]],],
              type="coefficients",
              s=lambda_choice
            ))[-1], n_responses, n_predictors, byrow=TRUE)
          })
        })

      }

    } else if (weightest == "multivar") {

      # Fit a standard lasso multivar model to get structured initial estimates
      # This is particularly useful for k=1 TVP where it properly separates
      # common from time-varying effects through joint structured estimation

      # Check if breaks are suitable for multivar weightest
      # Default breaks (single-timepoint periods) are not suitable
      if (tvp && length(breaks) > 0) {
        # Check if breaks appear to be the default (too many small periods)
        min_period_len <- min(sapply(breaks, function(b) min(sapply(b, length))))
        num_periods <- length(breaks[[1]])

        if (min_period_len < 10 || num_periods > 100) {
          stop(paste0(
            "weightest = 'multivar' requires explicit breaks with reasonable period sizes ",
            "for TVP models. Current breaks have ", num_periods, " periods with minimum ",
            "period length of ", min_period_len, ". Please provide explicit breaks via the ",
            "'breaks' argument, or use weightest = 'lasso' instead."
          ))
        }
      }

      # Build data list from Ak/bk
      # Reconstruct original time series: Ak[1,] is Y_1, bk is Y_2 to Y_T
      data_list <- lapply(seq_along(bk), function(i) {
        rbind(Ak[[i]][1, , drop=FALSE], bk[[i]])
      })

      # Convert breaks from window format back to break point format for constructModel
      # breaks[[i]] is a list of windows (index vectors), constructModel expects break points
      # Note: +1 adjustment because we prepended Ak[1,] to reconstruct the original data
      break_points_list <- lapply(seq_along(breaks), function(i) {
        windows <- breaks[[i]]
        if (length(windows) <= 1) {
          numeric(0)  # No breaks for single period
        } else {
          # Break points are at the end of each window except the last
          # Add 1 to account for the prepended row in reconstructed data
          sapply(windows[-length(windows)], function(w) max(w) + 1)
        }
      })

      # Construct a temporary model with standard lasso
      temp_model <- constructModel(
        data = data_list,
        tvp = tvp,
        breaks = break_points_list,
        lassotype = "standard",  # Standard lasso for initial estimates
        nlambda1 = 30,
        nfolds = nfolds,
        intercept = intercept,
        common_effects = common_effects,
        common_tvp_effects = common_tvp_effects
      )

      # Fit the model
      temp_fit <- cv.multivar(temp_model)

      # Extract the structured estimates
      total_effects <- temp_fit$mats$total

      # For TVP models, we get properly structured common and tvp estimates
      if (tvp) {
        common_effects_est <- temp_fit$mats$common
        unique_effects_est <- temp_fit$mats$unique
        tvp_effects_est <- temp_fit$mats$tvp
        common_tvp_effects_est <- temp_fit$mats$common_tvp

        # tvp_effects from breakup_transition is indexed by timepoint, not period
        # We need to extract period-level matrices for inittvpcoefs
        num_periods <- length(breaks[[1]])

        # Create inittvpcoefs: period-level total effects
        # This should be the FULL per-period estimate, not just TVP deviations
        inittvpcoefs <- lapply(seq_along(data_list), function(i) {
          lapply(1:num_periods, function(p) {
            # Get the first timepoint of period p to extract the period's TVP matrix
            first_t_in_period <- breaks[[i]][[p]][1]
            # Total per-period estimate = common + unique + tvp
            common_effects_est + unique_effects_est[[i]] + tvp_effects_est[[i]][[first_t_in_period]]
          })
        })

        # Compute tvp_effects as deviation from common
        # This ensures: common + tvp_effects[[i]][[p]] = inittvpcoefs[[i]][[p]]
        tvp_effects_computed <- lapply(seq_along(data_list), function(i) {
          lapply(1:num_periods, function(p) {
            inittvpcoefs[[i]][[p]] - common_effects_est
          })
        })

        # Return early for TVP case - skip the standard decomposition
        res <- list(
          common_effects = common_effects_est,
          subgroup_effects = NULL,
          unique_effects = unique_effects_est,
          tvp_effects = tvp_effects_computed,
          common_tvp_effects = common_tvp_effects_est,
          total_effects = total_effects,
          inittvpcoefs = inittvpcoefs
        )
        return(res)

      } else {
        # Non-TVP case: extract common and unique directly
        common_effects_est <- temp_fit$mats$common
        unique_effects_est <- temp_fit$mats$unique

        # Return early - skip the standard decomposition
        res <- list(
          common_effects = common_effects_est,
          subgroup_effects = NULL,
          unique_effects = unique_effects_est,
          tvp_effects = NULL,
          common_tvp_effects = NULL,
          total_effects = total_effects
        )
        return(res)
      }

    }

    # else if (weightest == "var") {
    # 
    #   if(tvp){
    #     
    #     inittvpcoefs <- lapply(seq_along(object@Ak),function(g){
    #       lapply(seq_along(breaks[[g]]), function(q){
    #         fit<- vars::VAR(as.matrix(object@Ak[[g]][breaks[[g]][[q]],]), p=1, type="none")$varresult
    #         as.matrix(do.call("rbind",lapply(seq_along(colnames(object@Ak[[g]])), function(x) {
    #           fit[[x]]$coefficients
    #         })))
    #       })
    #     })
    #     
    #   }
    #   
    # } 
    
    #---------------------------------------------#
    # 2. estimate common, unique, subgroup effects
    #---------------------------------------------#
    
    # if(weightest == "multivar2"){
    #   
    #   common_effects <- fit$mats$common
    #   unique_effects <- fit$mats$unique
    #   
    #   if(subgroup){
    #     
    #     subgroup_effects <- fit$mats$subgrp
    #     
    #   } else {
    #     
    #     subgroup_effects <- NULL
    #     
    #   }
    #   
    # } else {
    
      # create array of total effect matrices
      total_effects_array <- array(unlist(total_effects),
        c(dim(total_effects[[1]]), length(total_effects))
      )

      # Save the parameter value before overwriting
      include_common_effects <- common_effects

      if (tvp && !include_common_effects){
        # No common effects: skip median calculation
        common_effects <- NULL
      } else {
        # elementwise median of list of total effect matrices
        common_effects <- apply(total_effects_array, 1:2, median)
      }
      
      if(subgroup){
        
         subgroup_effects <- lapply(seq_along(1:max(subgroup_membership)), function(i){
           apply(total_effects_array[,,which(subgroup_membership==i)], 1:2, median)
         })
        
        # subgroup effects are deviation from common
        subgroup_effects <- lapply(seq_along(1:length(subgroup_effects)), function(i){
          subgroup_effects[[i]] <- subgroup_effects[[i]] - common_effects
        })

        # unique = total - common - subgroup ensures: total = common + subgroup + unique
        unique_effects <- lapply(seq_along(Ak), function(i){
          total_effects[[i]] - common_effects - subgroup_effects[[subgroup_membership[i]]]
        })
        
        tvp_effects <- NULL
        
      } else {
        
        subgroup_effects <- tvp_effects <-  NULL
        
        unique_effects <- lapply(seq_along(Ak), function(i){
          total_effects[[i]] - common_effects
        })
        
      }

      if (tvp) {

        subgroup_effects <- NULL

        # Special case: k=1 TVP
        # For k=1, there's no "common across subjects" - only periods
        if (k == 1) {

          num_periods <- length(inittvpcoefs[[1]])

          if (!include_common_effects) {
            # No time-invariant component: use per-period estimates directly
            common_effects <- NULL

            # unique_effects is zero for k=1
            unique_effects <- list(matrix(0, d, d))

            # TVP effects = per-period estimates directly (no common to subtract)
            tvp_effects <- list(inittvpcoefs[[1]])

          } else {
            # Derive common as median of per-period estimates
            # This selects edges that are consistently present across periods
            period_array <- array(unlist(inittvpcoefs[[1]]),
                                 c(dim(inittvpcoefs[[1]][[1]]), num_periods))
            common_effects <- apply(period_array, 1:2, median)

            # unique_effects is zero for k=1 (no subject-specific deviation)
            unique_effects <- list(matrix(0, nrow(common_effects), ncol(common_effects)))

            # TVP effects = period estimate - common
            # This ensures: common + tvp_effects[[p]] = inittvpcoefs[[p]]
            tvp_effects <- list(
              lapply(seq_along(inittvpcoefs[[1]]), function(p) {
                inittvpcoefs[[1]][[p]] - common_effects
              })
            )
          }

          common_tvp_effects_est <- NULL

        } else {
          # k > 1: original logic

          # Step 1: Compute common TVP effects (if enabled and k>1)
          if (common_tvp_effects && k > 1) {

            num_periods <- length(inittvpcoefs[[1]])

            # Pool estimates across subjects for each period
            common_tvp_effects_est <- lapply(1:num_periods, function(p) {
              # Collect period p estimates from all subjects
              period_estimates <- lapply(1:k, function(i) {
                inittvpcoefs[[i]][[p]]
              })

              # Take element-wise median across subjects
              period_array <- array(unlist(period_estimates),
                                   c(dim(period_estimates[[1]]), k))
              apply(period_array, 1:2, median)
            })

          } else {
            common_tvp_effects_est <- NULL
          }

          # Step 2: Compute base effects decomposition
          if (!include_common_effects){
            # No common effects: unique plays role of "base" (like common for k=1)
            # Keep unique = total for weights
            unique_effects <- total_effects

            # But for TVP, we want deviations from zero (like k=1)
            # Create a "pseudo-common" that's the average, to compute TVP correctly
            pseudo_common <- apply(total_effects_array, 1:2, median)
            pseudo_unique <- lapply(seq_along(Ak), function(i){
              total_effects[[i]] - pseudo_common
            })

          } else {
            # Standard logic: unique = total - common
            unique_effects <- lapply(seq_along(Ak), function(i){
              total_effects[[i]] - common_effects
            })
            pseudo_unique <- unique_effects  # For consistency below
          }

          # Step 3: Compute TVP effects decomposition
          if (common_tvp_effects && k > 1) {
            # Unique TVP = per-period estimate - common TVP
            tvp_effects <- lapply(seq_along(1:length(inittvpcoefs)), function(i){
              lapply(seq_along(1:length(inittvpcoefs[[i]])), function(j){
                inittvpcoefs[[i]][[j]] - common_tvp_effects_est[[j]]
              })
            })
          } else {
            # No common TVP: TVP = unique base - period estimate
            # (This is the current behavior)
            tvp_effects <- lapply(seq_along(1:length(inittvpcoefs)), function(i){
              lapply(seq_along(1:length(inittvpcoefs[[i]])), function(j){
                pseudo_unique[[i]] - inittvpcoefs[[i]][[j]]
              })
            })
          }

        }

      }
    
    #}
    
    # Add inittvpcoefs and common_tvp_effects to return list if tvp is TRUE
    if (tvp && exists("inittvpcoefs")) {
      res <- list(
        common_effects = common_effects,
        subgroup_effects = subgroup_effects,
        unique_effects = unique_effects,
        tvp_effects = tvp_effects,
        common_tvp_effects = if(exists("common_tvp_effects_est")) common_tvp_effects_est else NULL,
        total_effects = total_effects,
        inittvpcoefs = inittvpcoefs
      )
    } else {
      res <- list(
        common_effects = common_effects,
        subgroup_effects = subgroup_effects,
        unique_effects = unique_effects,
        tvp_effects = tvp_effects,
        common_tvp_effects = NULL,
        total_effects = total_effects
      )
    }
    
  }
  
  return(res)
  
}



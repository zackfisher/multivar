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
#' @param weightest Character. How to estimate initial coefficients for adaptive weights: "lasso", "ridge", or "ols". Only used when lassotype = "adaptive" (ignored for standard LASSO).
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
      
      
      # did not carry intercept through tvp yet
      if(tvp){

        inittvpcoefs <- lapply(seq_along(Ak),function(g){
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
            foldid_glmnet <- rep(period_foldid, ncol(bk[[g]]))

            lasso.fit <- glmnet::cv.glmnet(
              diag(ncol(Ak[[g]][breaks[[g]][[q]],])) %x% Ak[[g]][breaks[[g]][[q]],],
              as.vector(bk[[g]][breaks[[g]][[q]],]),
              family = "gaussian",
              alpha = 1,
              standardize = FALSE,
              foldid = foldid_glmnet
            )
            t(matrix(as.vector(predict(
              lasso.fit,
              diag(ncol(Ak[[g]][breaks[[g]][[q]],])) %x% Ak[[g]][breaks[[g]][[q]],],
              type="coefficients",
              s=lambda_choice
            ))[-1], ncol(Ak[[g]]), ncol(Ak[[g]])))
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



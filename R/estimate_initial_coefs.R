#' Estimate consistent first-stage coefficients
#' 
#' @description
#' `estimate_initial_coefs` returns consistent first-stage estimates to be used
#' for calculating adaptive weights used in regularization. 
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
    nfolds){
  
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
        
        g_coefs <- coef(lasso.fit, s = "lambda.1se")[-1]
        
        mat <- matrix(g_coefs, ncol(bk[[g]]), ncol(Ak[[g]]), byrow=TRUE)
        
        mat
        
      })
      
      
      # did not carry intercept through tvp yet
      if(tvp){

        inittvpcoefs <- lapply(seq_along(object@Ak),function(g){
          lapply(seq_along(breaks[[g]]), function(q){
            lasso.fit <- glmnet::cv.glmnet(
              diag(ncol(object@Ak[[g]][breaks[[g]][[q]],])) %x% object@Ak[[g]][breaks[[g]][[q]],], 
              as.vector(object@bk[[g]][breaks[[g]][[q]],]), 
              family = "gaussian", 
              alpha = 1, 
              standardize = FALSE, 
              nfolds = nfolds
            )
            t(matrix(as.vector(predict(
              lasso.fit, 
              diag(ncol(object@Ak[[g]][breaks[[g]][[q]],])) %x% object@Ak[[g]][breaks[[g]][[q]],], 
              type="coefficients", 
              s="lambda.1se"
            ))[-1], ncol(object@Ak[[g]]), ncol(object@Ak[[g]])))
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
        
        g_coefs <- coef(ridge.fit, s = "lambda.1se")[-1]
        
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
      
      # elementwise median of list of total effect matrices
      common_effects <- apply(total_effects_array, 1:2, median)
      
      if(subgroup){
        
         subgroup_effects <- lapply(seq_along(1:max(subgroup_membership)), function(i){
           apply(total_effects_array[,,which(subgroup_membership==i)], 1:2, median)
         })
        
        # zf: is this backwards?
        subgroup_effects <- lapply(seq_along(1:length(subgroup_effects)), function(i){
          subgroup_effects[[i]] <- subgroup_effects[[i]] - common_effects
        })
        
        unique_effects <- lapply(seq_along(Ak), function(i){
          # already subtracted out the common effects
          # total_effects[[i]] - common_effects - subgroup_effects[[subgroup[i]]]
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
        
        unique_effects <- lapply(seq_along(Ak), function(i){
          total_effects[[i]] - common_effects
        })
        
        tvp_effects <- lapply(seq_along(1:length(inittvpcoefs)), function(i){
          lapply(seq_along(1:length(inittvpcoefs[[i]])), function(j){
            unique_effects[[i]] - inittvpcoefs[[i]][[j]]
          })
        })
        
        
      }
    
    #}
    
    res <- list(
      common_effects = common_effects,
      subgroup_effects = subgroup_effects,
      unique_effects = unique_effects,
      tvp_effects = tvp_effects,
      total_effects = total_effects
    )
    
  }
  
  return(res)
  
}



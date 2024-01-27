#' Estimate consistent first-stage coefficients
#' 
#' @description
#' `estimate_initial_coefs` returns consistent first-stage estimates to be used
#' for calculating adaptive weights used in regularization. 
#'
#' @details
#' This function returns a list containing the following elements:
#' * common_effects: a d x d matrix containing consistent esimtates of the common 
#'   effects
#' * subgroup_effects: a list of length max(subgroup), where each element is a d x d matrix 
#'   containing consistent esimtates of the subgroup effects
#' * unique_effects: a list of length k, where each element is a d x d matrix 
#'   containing consistent esimtates of the unique effects
#' * total_effects: a list of length k, where each element is a d x d matrix 
#'   containing consistent esimtates of the total effects
#' 
#' #' @export
estimate_initial_coefs <- function(
    Ak, 
    bk,
    d, 
    k, 
    lassotype, 
    weightest,
    subgroup,
    subgroupflag,
    nlambda1,
    nlambda2){
  
  # Ak<-object@Ak
  # bk<- object@bk
  # ratios<-object@ratios
  # d<-object@d
  # k<-object@k
  # lassotype<-object@lassotype
  # weightest<-object@weightest
  # subgroup <- object@subgroup
  # subgroupflag <- object@subgroupflag
  # nlambda1 <- object@nlambda1
  # nlambda2 <- object@nlambda2
  
  if (lassotype == "standard"){
    
    if(!subgroupflag){
      res <- list(
        common_effects = matrix(1, d, d),
        subgroup_effects = NULL,
        unique_effects = replicate(k,  matrix(1, d, d)),
        total_effects = replicate(k,  matrix(1, d, d))
      )
    } else {
      res <- list(
        common_effects = matrix(1, d, d),
        subgroup_effects = replicate(max(subgroup),  matrix(1, d, d)),
        unique_effects = replicate(k,  matrix(1, d, d)),
        total_effects = replicate(k,  matrix(1, d, d))
      )
    }
    
  } else {
    
    #---------------------------------------------#
    # 1. estimate total effects
    #---------------------------------------------#
    
    if (weightest == "ols") {
        
      total_effects <- lapply(seq_along(Ak), function(z){
        b_ols_i <- matrix(0, ncol(Ak[[z]]),ncol(Ak[[z]]))
        for(i in 1:ncol(Ak[[z]])){
          for(j in 1:ncol(Ak[[z]])){
            b_ols_i[j,i] <- lm.fit(as.matrix(Ak[[z]])[,i,drop=F],as.matrix(bk[[z]])[,j,drop=F])$coefficients
          }
        }
        b_ols_i
      })
        
    } else if (weightest == "lasso") {
      
      total_effects <- lapply(seq_along(Ak),function(g){
        lasso.fit <- glmnet::cv.glmnet(
          diag(ncol(Ak[[g]])) %x% Ak[[g]], 
          as.vector(bk[[g]]), 
          family = "gaussian", 
          alpha = 1, 
          standardize = FALSE, 
          nfolds = 10
        )
        t(matrix(as.vector(predict(
          lasso.fit, 
          diag(ncol(Ak[[g]])) %x% Ak[[g]], 
          type="coefficients", 
          s="lambda.1se"
        ))[-1], ncol(Ak[[g]]), ncol(Ak[[g]])))
      })
      
    } else if (weightest == "ridge") {
      
      total_effects <- lapply(seq_along(Ak),function(g){
        lasso.fit <- glmnet::cv.glmnet(
          diag(ncol(Ak[[g]])) %x% Ak[[g]], 
          as.vector(bk[[g]]), 
          family = "gaussian", 
          alpha = 0, 
          standardize = FALSE, 
          nfolds = 10
        )
        t(matrix(as.vector(predict(
          lasso.fit, 
          diag(ncol(Ak[[g]])) %x% Ak[[g]], 
          type="coefficients", 
          s="lambda.1se"
        ))[-1], ncol(Ak[[g]]), ncol(Ak[[g]])))
      })
      
    } else if (weightest == "var") {
      
      total_effects <- lapply(seq_along(Ak),function(g){
        fit<- vars::VAR(as.matrix(Ak[[g]]), p=1, type="none")$varresult
        as.matrix(do.call("rbind",lapply(seq_along(colnames(Ak[[g]])), function(x) {
          fit[[x]]$coefficients
        })))
      })
      
    } else if (weightest == "multivar1" | weightest == "multivar2" ) {
      
      if(subgroupflag){
        
        mod <- constructModel(
          data = Ak, 
          lassotype = "standard",
          nlambda1 = nlambda1, 
          nlambda2 = nlambda2,
          subgroup = subgroup)
        
        fit <- cv.multivar(mod)
        total_effects <- fit$mats$total
        
      }else{
        
        mod <- constructModel(
          data = Ak, 
          lassotype = "standard",
          nlambda1 = nlambda1, 
          nlambda2 = nlambda2
        )
        
        fit <- cv.multivar(mod)
        total_effects <- fit$mats$total
        
      }
 
    } else if (weightest == "variable") {
      
      if(ncol(Ak[[1]]) >= nrow(Ak[[1]])){
        
        total_effects <- lapply(seq_along(Ak), function(z){
          b_ols_i <- matrix(0, ncol(Ak[[z]]),ncol(Ak[[z]]))
          for(i in 1:ncol(Ak[[z]])){
            for(j in 1:ncol(Ak[[z]])){
              b_ols_i[j,i] <- lm.fit(as.matrix(Ak[[z]])[,i,drop=F],as.matrix(bk[[z]])[,j,drop=F])$coefficients
            }
          }
          b_ols_i
        })
        
      } else {
        
        total_effects <- lapply(seq_along(Ak),function(g){
          fit<- vars::VAR(as.matrix(Ak[[g]]), p=1, type="none")$varresult
          as.matrix(do.call("rbind",lapply(seq_along(colnames(Ak[[g]])), function(x) {
            fit[[x]]$coefficients
          })))
        })
        
      }
      
    }
    
    #---------------------------------------------#
    # 2. estimate common, unique, subgroup effects
    #---------------------------------------------#
    
    if(weightest == "multivar2"){
      
      common_effects <- fit$mats$common
      unique_effects <- fit$mats$unique
      
      if(subgroupflag){
        
        subgroup_effects <- fit$mats$subgroup
        
      } else {
        
        subgroup_effects <- NULL
        
      }
      
    } else {
    
      # create array of total effect matrices
      total_effects_array <- array(unlist(total_effects), 
        c(dim(total_effects[[1]]), length(total_effects))
      )
      
      # elementwise median of list of total effect matrices
      common_effects <- apply(total_effects_array, 1:2, median)
      
      if(subgroupflag){
        
        subgroup_effects <- lapply(seq_along(1:max(subgroup)), function(i){
          apply(total_effects_array[,,which(subgroup==i)], 1:2, median)
        })
        
        # zf: is this backwards?
        subgroup_effects <- lapply(seq_along(1:length(subgroup_effects)), function(i){
          subgroup_effects[[i]] <- subgroup_effects[[i]] - common_effects
        })
        
        unique_effects <- lapply(seq_along(Ak), function(i){
          # already subtracted out the common effects
          # total_effects[[i]] - common_effects - subgroup_effects[[subgroup[i]]]
          total_effects[[i]] - subgroup_effects[[subgroup[i]]]
        })
        
      } else {
        
        subgroup_effects <- NULL
        
        unique_effects <- lapply(seq_along(Ak), function(i){
          total_effects[[i]] - common_effects
        })
        
      }
    
    }
    
    res <- list(
      common_effects = common_effects,
      subgroup_effects = subgroup_effects,
      unique_effects = unique_effects,
      total_effects = total_effects
    )
    
  }
  
  return(res)
  
}



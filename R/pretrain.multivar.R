#' @export
pretrain.multivar <- function(object, alpha = 0, lambda_best = "1se", stage1adaptive = FALSE, stage2adaptive = FALSE) {
  

  stage1 <- pretrained_var_stage1(
    object = object, 
    stage1adaptive = stage1adaptive, 
    stage1ratios = FALSE, 
    lambda_best = lambda_best
  )
  
  #stage1_glmnet <- pretrained_var_stage1_glmnet(object)

  stage2 <- pretrained_var_stage2(
    object = object, 
    stage1 = stage1, 
    alpha = alpha, 
    stage2adaptive =stage2adaptive, 
    stage2ratios = FALSE, 
    intercept = FALSE, 
    lambda_best = lambda_best
  ) 
      
  if(object@intercept){
    common <- stage1$A_pre[,-1,drop=F]
    unique <- lapply(stage2,function(x){ x$U })
    total  <- lapply(unique$A,function(x){ x$A })
  } else {
    common <- stage1$A_pre[,-1,drop=F]
    unique <- lapply(stage2,function(x){ x$U[,-1,drop=F]})
    total  <- lapply(stage2,function(x){ x$A[,-1,drop=F]})
  }
 
  
  mats <- list(
   "common" = common,
   "subgrp" = NULL,
   "unique" = unique,
   "total" = total
  )
  
  return(list(
    stage1 = stage1,
    stage2 = stage2,
    MSFE_min = c(stage1$MSFE_min,lapply(stage2,function(x){x$MSFE_min})),
    mats = mats
  ))

}
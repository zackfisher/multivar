#' @keywords internal
cv.pretrain.multivar <- function(object, lambda_best = "min", stage1adaptive = TRUE, stage2adaptive = TRUE) {
  
  alpha <- seq(0,1,by=0.1)
  res <- list()
  for(i in 1:length(alpha)){
    res[[i]] <- pretrain.multivar(
      object, 
      alpha = alpha[i],
      lambda_best = lambda_best, 
      stage1adaptive = stage1adaptive, 
      stage2adaptive = stage2adaptive
    )
  }
  
  optimal_alpha_idx <- which.min(lapply(res,function(x){
    mean( unlist(x$MSFE_min))
  }))
  
  #optimal_alpha_idx <- 6
  
  return(list(
    full = res,
    stage1 = res[[optimal_alpha_idx]]$stage1,
    stage2 = res[[optimal_alpha_idx]]$stage2,
    #MSFE = stage2$MSFE,
    mats = res[[optimal_alpha_idx]]$mats
  ))

}
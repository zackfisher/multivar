#' Extract hyperparameter info from multivar CV fit
#'
#' Run this right after:
#'   fit <- cv_multivar(...)
#'
#' It reconstructs:
#'   - full lambda1 grid (from the model object)
#'   - implied lambda2 grid (lambda1 * ratio)
#'   - scenario table (one row per scenario)
#'   - MSFE surface (nlam x n_scenarios) from fit[[2]]
#'   - best hyperparameter combo (lambda1, ratio, lambda2)
#'
#' @param object multivar model object used to call cv_multivar(...)
#' @param fit    result returned by cv_multivar(...); MSFE is assumed to be fit[[2]]
#' @export
extract_multivar_hyperparams <- function(object, fit) {
  
  
  
  
  # number of lambda1s and ratios
  n_lambda <- nrow(object@lambda1)
  n_ratios <- ncol(object@lambda1)
  
  # Create all combinations of lambda1 and ratio indices
  param_map <- expand.grid(
    ratio_index = seq_len(n_lambda),
    lambda1_index   = seq_len(n_ratios)
  )
  
  # Compute slice index consistent with C++ ordering: i*n_r + j
  param_map$slice_index <- (param_map$lambda1_index - 1) * n_ratios + (param_map$ratio_index - 1)
  
  # Assign the actual lambda1 and ratio values
  param_map$lambda1_value <- as.vector(t(object@lambda1))
  param_map$ratio_value   <- object@ratios[param_map$ratio_index]
  
  
  MSFE_vec <- as.vector(t(colMeans(fit[[2]])))  # ensures order matches i*n_r + j in C++
  
  cv_results <- param_map
  cv_results$MSFE <- MSFE_vec
  
  return(cv_results)
  
  ## ------------------------------------------------------------
  ## 1. lambda1 grid (must already be in the model object)
  ## ------------------------------------------------------------
  # lambda1_grid <- object@lambda1
  # if (is.null(lambda1_grid)) {
  #   stop("extract_multivar_hyperparams: object@lambda1 is NULL.")
  # }
  # nlam  <- nrow(lambda1_grid)
  # nscen <- ncol(lambda1_grid)
  # 
  # ## ------------------------------------------------------------
  # ## 2. scenario table (you already added build_scenario_table)
  # ## ------------------------------------------------------------
  # scen_tab <- build_scenario_table(object)
  # if (nrow(scen_tab) != nscen) {
  #   warning("extract_multivar_hyperparams: scenario_table rows (",
  #           nrow(scen_tab),
  #           ") != ncol(lambda1_grid) (",
  #           nscen, "). Check scenario ordering.")
  # }
  # 
  # ## ------------------------------------------------------------
  # ## 3. lambda2 grid = lambda1 * ratio_scenario
  # ##    (simple, no-subgroup/no-tvp case)
  # ## ------------------------------------------------------------
  # lambda2_grid <- sweep(lambda1_grid, 2, scen_tab$ratio, "*")
  # 
  # ## ------------------------------------------------------------
  # ## 4. get MSFE from fit[[2]]
  # ##    your current cv_multivar returns CV errors in position 2
  # ## ------------------------------------------------------------
  # msfe_raw <- fit[[2]]
  # if (is.null(msfe_raw)) {
  #   stop("extract_multivar_hyperparams: fit[[2]] is NULL; can't find MSFE.")
  # }
  # 
  # ## ------------------------------------------------------------
  # ## 5. reshape MSFE to (nlam x nscen)
  # ##    common pattern: msfe_raw is folds x (nlam * nscen)
  # ## ------------------------------------------------------------
  # if (is.matrix(msfe_raw) && ncol(msfe_raw) == nlam * nscen) {
  #   # case: folds x combos
  #   msfe_vec  <- colMeans(msfe_raw)
  #   msfe_grid <- matrix(msfe_vec, nrow = nlam, ncol = nscen)
  # } else if (is.vector(msfe_raw) && length(msfe_raw) == nlam * nscen) {
  #   # case: flat vector
  #   msfe_grid <- matrix(msfe_raw, nrow = nlam, ncol = nscen)
  # } else if (is.matrix(msfe_raw) &&
  #            all(dim(msfe_raw) == c(nlam, nscen))) {
  #   # already the right shape
  #   msfe_grid <- msfe_raw
  # } else {
  #   stop("extract_multivar_hyperparams: can't reshape fit[[2]] into (nlam x n_scen). ",
  #        "Got dim = ", paste(dim(msfe_raw), collapse = " x "),
  #        ", expected ", nlam, " x ", nscen, " (possibly with folds).")
  # }
  # 
  # ## ------------------------------------------------------------
  # ## 6. find best (lambda1 row, scenario column)
  # ## ------------------------------------------------------------
  # best_idx <- which(msfe_grid == min(msfe_grid, na.rm = TRUE),
  #                   arr.ind = TRUE)[1, ]
  # best_row <- best_idx[1]       # which lambda1 along the path
  # best_col <- best_idx[2]       # which scenario
  # 
  # best_lambda1 <- lambda1_grid[best_row, best_col]
  # best_ratio   <- scen_tab$ratio[best_col]
  # best_lambda2 <- best_lambda1 * best_ratio
  # 
  # ## ------------------------------------------------------------
  # ## 7. return everything
  # ## ------------------------------------------------------------
  # list(
  #   lambda1_grid   = lambda1_grid,
  #   lambda2_grid   = lambda2_grid,
  #   scenario_table = scen_tab,
  #   msfe_grid      = msfe_grid,
  #   msfe_raw       = msfe_raw,
  #   best_row       = best_row,
  #   best_scenario  = best_col,
  #   best_lambda1   = best_lambda1,
  #   best_ratio     = best_ratio,
  #   best_lambda2   = best_lambda2
  # )
}

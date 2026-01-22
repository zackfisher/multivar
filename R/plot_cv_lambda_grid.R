#' Plot CV error over (lambda1, lambda2)
#'
#' @param fit    result from cv.multivar(...)
#' @param object OPTIONAL; only used if fit does NOT contain the original model as fit$obj
#' @param log_scale logical; plot on log10 scale?
#' @param show_best logical; mark best pair?
#' @export
plot_cv_lambda_grid <- function(fit,
                                object = NULL,
                                log_scale = TRUE,
                                show_best = TRUE) {
  
  # 1. try to use hyperparams already stored in fit
  if (!is.null(fit$hyperparams)) {
    hyp <- fit$hyperparams
    
  } else {
    # 2. try to get the original model from fit$obj
    if (is.null(object)) {
      if (!is.null(fit$obj)) {
        object <- fit$obj
      } else {
        stop("plot_cv_lambda_grid: neither fit$hyperparams nor fit$obj was found, ",
             "and you didn't supply 'object'.")
      }
    }
    
    # 3. fall back to your existing extractor
    hyp <- extract_multivar_hyperparams(object, fit)
  }
  
  lambda1_grid <- hyp$lambda1_grid
  lambda2_grid <- hyp$lambda2_grid
  msfe_grid    <- hyp$msfe_grid
  
  nlam  <- nrow(lambda1_grid)
  nscen <- ncol(lambda1_grid)
  
  df <- data.frame(
    lambda1 = as.numeric(lambda1_grid),
    lambda2 = as.numeric(lambda2_grid),
    msfe    = as.numeric(msfe_grid),
    row     = rep(seq_len(nlam), times = nscen),
    scen    = rep(seq_len(nscen), each  = nlam)
  )
  
  if (log_scale) {
    df$log_lambda1 <- log10(df$lambda1)
    df$log_lambda2 <- log10(df$lambda2)
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not installed; returning grid data.")
    return(invisible(df))
  }
  library(ggplot2)
  
  if (log_scale) {
    p <- ggplot(df, aes(x = log_lambda1, y = log_lambda2, fill = msfe)) +
      geom_tile() +
      scale_fill_viridis_c(option = "magma", direction = -1) +
      labs(
        x = "log10(lambda1)",
        y = "log10(lambda2)",
        fill = "CV error",
        title = "CV error over (lambda1, lambda2)"
      ) +
      theme_minimal()
  } else {
    p <- ggplot(df, aes(x = lambda1, y = lambda2, fill = msfe)) +
      geom_tile() +
      scale_fill_viridis_c(option = "magma", direction = -1) +
      labs(
        x = "lambda1",
        y = "lambda2",
        fill = "CV error",
        title = "CV error over (lambda1, lambda2)"
      ) +
      theme_minimal()
  }
  
  if (show_best && !is.null(hyp$best_row) && !is.null(hyp$best_scenario)) {
    best_df <- df[df$row == hyp$best_row & df$scen == hyp$best_scenario, , drop = FALSE]
    if (nrow(best_df) == 1) {
      if (log_scale) {
        p <- p + geom_point(
          data  = best_df,
          aes(x = log_lambda1, y = log_lambda2),
          size  = 3,
          shape = 21,
          fill  = "black",
          color = "cyan"
        )
      } else {
        p <- p + geom_point(
          data  = best_df,
          aes(x = lambda1, y = lambda2),
          size  = 3,
          shape = 21,
          fill  = "black",
          color = "cyan"
        )
      }
    }
  }
  
  print(p)
  invisible(df)
}

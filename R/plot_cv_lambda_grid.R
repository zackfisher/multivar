#' Plot CV error over (lambda1, lambda2)
#'
#' @param fit    result from cv.multivar(...)
#' @param object OPTIONAL; only used if fit does NOT contain hyperparams or obj
#' @param log_scale logical; plot on log10 scale?
#' @param show_best logical; mark best pair?
#' @export
plot_cv_lambda_grid <- function(fit,
                                object = NULL,
                                log_scale = TRUE,
                                show_best = TRUE) {
  
  # 1) Get hyperparams dataframe
  if (!is.null(fit$hyperparams)) {
    df <- fit$hyperparams
  } else {
    if (is.null(object)) {
      if (!is.null(fit$obj)) {
        object <- fit$obj
      } else {
        stop("plot_cv_lambda_grid: neither fit$hyperparams nor fit$obj was found, and you didn't supply 'object'.")
      }
    }
    df <- extract_multivar_hyperparams(object, fit)
  }
  
  # 2) Prepare plotting dataframe
  if (log_scale) {
    df$log_lambda1 <- log10(df$lambda1_value)
    df$log_lambda2 <- log10(df$lambda2_value)
  }
  
  # 3) Find best hyperparameters
  best_idx <- which.min(df$MSFE)
  best_row <- df[best_idx, , drop = FALSE]
  
  # 4) Check for ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not installed; returning grid data.")
    return(invisible(df))
  }
  
  # Optional viridis scale (fallback if not installed)
  has_viridis <- requireNamespace("viridis", quietly = TRUE)
  
  # 5) Create plot (use .data pronoun to satisfy R CMD check)
  if (log_scale) {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$log_lambda1,
        y = .data$log_lambda2,
        fill = .data$MSFE
      )
    ) +
      ggplot2::geom_tile() +
      (if (has_viridis) {
        viridis::scale_fill_viridis(discrete = FALSE, option = "magma", direction = -1)
      } else {
        ggplot2::scale_fill_gradient()
      }) +
      ggplot2::labs(
        x = "log10(lambda1)",
        y = "log10(lambda2)",
        fill = "CV error",
        title = "CV error over (lambda1, lambda2)"
      ) +
      ggplot2::theme_minimal()
  } else {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$lambda1_value,
        y = .data$lambda2_value,
        fill = .data$MSFE
      )
    ) +
      ggplot2::geom_tile() +
      (if (has_viridis) {
        viridis::scale_fill_viridis(discrete = FALSE, option = "magma", direction = -1)
      } else {
        ggplot2::scale_fill_gradient()
      }) +
      ggplot2::labs(
        x = "lambda1",
        y = "lambda2",
        fill = "CV error",
        title = "CV error over (lambda1, lambda2)"
      ) +
      ggplot2::theme_minimal()
  }
  
  # 6) Mark best point
  if (show_best && nrow(best_row) == 1) {
    if (log_scale) {
      p <- p + ggplot2::geom_point(
        data = best_row,
        ggplot2::aes(
          x = .data$log_lambda1,
          y = .data$log_lambda2
        ),
        size = 3,
        shape = 21,
        fill = "black",
        color = "cyan"
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = best_row,
        ggplot2::aes(
          x = .data$lambda1_value,
          y = .data$lambda2_value
        ),
        size = 3,
        shape = 21,
        fill = "black",
        color = "cyan"
      )
    }
  }
  
  print(p)
  invisible(df)
}

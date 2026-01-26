#' Plot fitted multivar model results
#'
#' Creates heatmap visualizations of transition matrices from fitted multivar models.
#' Supports plotting common effects, unique (subject-specific) effects, total effects,
#' subgroup effects, and time-varying parameter (TVP) effects.
#'
#' @param x Object returned by \code{cv.multivar()} or related fitting functions
#' @param plot_type Character. Type of effects to plot:
#'   \itemize{
#'     \item "common" - Common effects shared across all subjects
#'     \item "unique" - Subject-specific unique effects
#'     \item "total" - Total effects (common + unique for each subject)
#'     \item "subgrp" - Subgroup effects (if model includes subgroups)
#'     \item "tvp" - Time-varying parameter effects (unique TVP deviations)
#'     \item "common_tvp" - Common time-varying effects across subjects
#'     \item "tvp_total" - Total effects for TVP models by period
#'   }
#' @param facet_ncol Numeric. Number of columns when faceting multiple matrices
#' @param subjects Character "all" or numeric vector. Which subjects to plot for
#'   "unique", "total", or "tvp" types. Default is "all"
#' @param periods Character "all" or numeric vector. Which time periods to plot for
#'   TVP models. Default is "all"
#' @param ub Numeric. Upper bound for color scale. Default is 1
#' @param lb Numeric. Lower bound for color scale. Default is -1
#' @param palette Character. Color palette: "default", "viridis", or "greyscale"
#' @param show_zeros Logical. If FALSE (default), zero values shown as white/NA
#' @param ... Additional arguments passed to plotting engine
#'
#' @return A ggplot2 object that can be further customized
#'
#' @examples
#' \dontrun{
#' sim <- multivar_sim(k = 3, d = 5, n = 50, prop_fill_com = 0.2,
#'                     prop_fill_ind = 0.1, sigma = diag(5))
#' model <- constructModel(sim$data)
#' fit <- cv.multivar(model)
#'
#' # Plot common effects
#' plot_results(fit, plot_type = "common")
#'
#' # Plot unique effects for first 2 subjects
#' plot_results(fit, plot_type = "unique", subjects = 1:2)
#'
#' # For TVP models
#' fit_tvp <- cv.multivar(constructModel(data, tvp = TRUE, breaks = list(c(50, 100))))
#' plot_results(fit_tvp, plot_type = "tvp")  # TVP deviations by period
#' plot_results(fit_tvp, plot_type = "tvp_total")  # Total effects by period
#'
#' # Customize the plot
#' p <- plot_results(fit, plot_type = "total")
#' p + ggplot2::ggtitle("My Custom Title") + ggplot2::theme_minimal()
#' }
#'
#' @seealso \code{\link{plot_sim}}, \code{\link{plot_transition_mat}}
#' @export
plot_results <- function(x,
                        plot_type = c("common", "unique", "total", "subgrp", "tvp", "common_tvp", "tvp_total"),
                        facet_ncol = 3,
                        subjects = "all",
                        periods = "all",
                        ub = 1,
                        lb = -1,
                        palette = "default",
                        show_zeros = FALSE,
                        ...) {

  # Validate plot_type
  plot_type <- match.arg(plot_type)

  # Handle tvp_total as a special case (plot total effects by period)
  if (plot_type == "tvp_total") {
    # For TVP models, total effects are already organized by period
    # For K=1: mats$total[[1]] is a list of period matrices
    # For K>1: mats$total[[i]] is a list of time-point matrices

    mats <- x$mats
    if (is.null(mats$total)) {
      stop("No total effects in this model")
    }

    k <- length(mats$total)

    if (k == 1) {
      # K=1: Total effects organized by period
      period_mats <- mats$total[[1]]

      # Select periods
      if (periods[1] == "all") {
        period_idx <- seq_along(period_mats)
      } else {
        period_idx <- periods
      }

      mat_list <- period_mats[period_idx]
      names(mat_list) <- paste0("Total Period ", period_idx)

    } else {
      # K>1: Select subjects first
      if (subjects[1] == "all") {
        subject_idx <- seq_along(mats$total)
      } else {
        subject_idx <- subjects
      }

      # For each subject, extract period-specific matrices
      mat_list <- list()
      for (s in subject_idx) {
        if (is.list(mats$total[[s]])) {
          num_periods <- length(mats$total[[s]])

          # Select periods
          if (periods[1] == "all") {
            period_idx <- seq_len(num_periods)
          } else {
            period_idx <- periods
          }

          for (p in period_idx) {
            mat_list[[length(mat_list) + 1]] <- mats$total[[s]][[p]]
            names(mat_list)[length(mat_list)] <- paste0("Subject ", s, " Period ", p)
          }
        }
      }
    }

    title <- "Total Transition Matrices by Period"

  } else {
    # Standard extraction for non-TVP types
    mat_list <- .extract_matrices_from_fit(
      fit = x,
      type = plot_type,
      subjects = subjects,
      periods = periods
    )

    # Determine title
    title <- switch(plot_type,
      common = "Common Transition Matrix",
      unique = "Unique Transition Matrices",
      total = "Total Transition Matrices",
      subgrp = "Subgroup Transition Matrices",
      tvp = "TVP Transition Matrices",
      common_tvp = "Common TVP Transition Matrices"
    )
  }

  # Call core plotting engine
  .plot_transition_heatmap(
    mat_list = mat_list,
    titles = names(mat_list),
    facet_ncol = facet_ncol,
    lb = lb,
    ub = ub,
    show_zeros = show_zeros,
    palette = palette,
    title = title,
    ...
  )
}

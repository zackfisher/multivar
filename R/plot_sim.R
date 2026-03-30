#' Plot simulated multivar ground truth
#'
#' Creates heatmap visualizations of true transition matrices from multivar simulations.
#' Supports plotting common effects, unique (subject-specific) effects, total effects,
#' and subgroup effects.
#'
#' @param x Object returned by \code{multivar_sim()} or \code{multivar_sim_subgroups()}
#' @param plot_type Character. Type of effects to plot:
#'   \itemize{
#'     \item "common" - Common effects shared across all subjects
#'     \item "unique" - Subject-specific unique effects
#'     \item "total" - Total effects (common + unique for each subject)
#'     \item "subgrp" - Subgroup effects (if simulation includes subgroups)
#'   }
#' @param facet_ncol Numeric. Number of columns when faceting multiple matrices
#' @param subjects Character "all" or numeric vector. Which subjects to plot for
#'   "unique" or "total" types. Default is "all"
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
#' # Simulate data
#' sim <- multivar_sim(k = 3, d = 5, n = 50,
#'                     prop_fill_com = 0.2, prop_fill_ind = 0.1,
#'                     lb = 0.1, ub = 0.9, sigma = diag(5))
#'
#' # Plot true common effects
#' plot_sim(sim, plot_type = "common")
#'
#' # Plot true unique effects for subjects 1-2
#' plot_sim(sim, plot_type = "unique", subjects = 1:2)
#'
#' # Simulate with subgroups
#' sim_sub <- multivar_sim_subgroups(k = 4, d = 3, n = 40,
#'                                   subgroup = c(1, 1, 2, 2),
#'                                   p_com = 0.25, p_sub = 0.10, p_ind = 0.05,
#'                                   sigma = diag(3))
#'
#' # Plot subgroup effects
#' plot_sim(sim_sub, plot_type = "subgrp")
#' }
#'
#' @seealso \code{\link{plot_results}}, \code{\link{plot_transition_mat}}
#' @export
plot_sim <- function(x,
                    plot_type = c("common", "unique", "total", "subgrp"),
                    facet_ncol = 3,
                    subjects = "all",
                    ub = 1,
                    lb = -1,
                    palette = "default",
                    show_zeros = FALSE,
                    ...) {

  # Validate plot_type
  plot_type <- match.arg(plot_type)

  # Extract matrices using helper function
  mat_list <- .extract_matrices_from_sim(
    sim = x,
    type = plot_type,
    subjects = subjects
  )

  # Determine title
  title <- switch(plot_type,
    common = "True Common Transition Matrix",
    unique = "True Unique Transition Matrices",
    total = "True Total Transition Matrices",
    subgrp = "True Subgroup Transition Matrices"
  )

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

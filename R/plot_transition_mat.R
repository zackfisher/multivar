#' Plot arbitrary transition matrix
#'
#' Creates a heatmap visualization of an arbitrary transition matrix.
#' Useful for plotting custom matrices or comparing different estimates.
#'
#' @param x Matrix. A transition matrix to plot
#' @param title Character. Main title for the plot
#' @param subtitle Character. Subtitle for the plot
#' @param ub Numeric. Upper bound for color scale. Default is 1
#' @param lb Numeric. Lower bound for color scale. Default is -1
#' @param legend Logical. Should a legend be included? Default TRUE
#' @param dimnames Logical. Should variable names be shown? Default TRUE
#' @param palette Character. Color palette: "default", "viridis", or "greyscale"
#' @param show_zeros Logical. If FALSE (default), zero values shown as white/NA
#' @param ... Additional arguments passed to plotting engine
#'
#' @return A ggplot2 object that can be further customized
#'
#' @examples
#' \dontrun{
#' # Plot random matrix
#' plot_transition_mat(matrix(rnorm(25), 5, 5), title = "Random Matrix")
#'
#' # Plot with custom bounds
#' mat <- matrix(runif(25, -0.5, 0.5), 5, 5)
#' plot_transition_mat(mat, title = "Small Coefficients", lb = -0.5, ub = 0.5)
#'
#' # Without legend or dimension names
#' plot_transition_mat(mat, legend = FALSE, dimnames = FALSE)
#' }
#'
#' @seealso \code{\link{plot_results}}, \code{\link{plot_sim}}
#' @export
plot_transition_mat <- function(x,
                               title = NULL,
                               subtitle = NULL,
                               ub = 1,
                               lb = -1,
                               legend = TRUE,
                               dimnames = TRUE,
                               palette = "default",
                               show_zeros = FALSE,
                               ...) {

  if (!is.matrix(x)) {
    stop("x must be a matrix")
  }

  # Remove dimnames if requested
  if (!dimnames) {
    colnames(x) <- NULL
    rownames(x) <- NULL
  }

  # Call core plotting engine
  p <- .plot_transition_heatmap(
    mat_list = x,
    titles = NULL,
    facet_ncol = 1,
    lb = lb,
    ub = ub,
    show_zeros = show_zeros,
    palette = palette,
    show_dimnames = dimnames,
    title = title,
    subtitle = subtitle,
    ...
  )

  # Remove legend if requested
  if (!legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  return(p)
}

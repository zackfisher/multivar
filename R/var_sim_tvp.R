#' Simulate a time-varying VAR model (piecewise or continuous)
#'
#' @param n Integer. Total time series length.
#' @param sigma Matrix. Innovation covariance matrix (d x d).
#' @param phi_list List of d x d transition matrices. For piecewise: one per period.
#'   For continuous: one base matrix.
#' @param T_per_period Integer vector. Observations per period (piecewise only).
#'   If NULL, periods are equal length.
#' @param mat_tvp Logical matrix (d x d). Which elements are time-varying (continuous only).
#' @param type Character. Either "piecewise" or "continuous".
#' @param burn Integer. Burn-in period to discard. Default 500.
#' @param intercept Intercept term. Can be:
#'   - Scalar (same for all variables and periods)
#'   - Vector of length d (different per variable, constant over periods)
#'   - List of vectors (different per period, for piecewise TVP)
#'
#' @return List with components:
#'   \item{data}{List containing the T x d data matrix}
#'   \item{mat_ind_final}{List containing the true transition matrices (list of period matrices)}
#'   \item{T_per_period}{Observations per period}
#'   \item{breaks}{Break points for period changes}
#'   \item{intercept}{Intercept values used}
#'
#' @examples
#' # Piecewise TVP
#' phi1 <- matrix(c(0.3, 0, 0, 0.3), 2, 2)
#' phi2 <- matrix(c(0.3, 0.2, 0, 0.3), 2, 2)
#' sim <- var_sim_tvp(n = 200, sigma = diag(2), phi_list = list(phi1, phi2),
#'                    T_per_period = c(100, 100), type = "piecewise")
#'
#' # Piecewise with period-specific intercepts
#' sim <- var_sim_tvp(n = 200, sigma = diag(2), phi_list = list(phi1, phi2),
#'                    T_per_period = c(100, 100), type = "piecewise",
#'                    intercept = list(c(1, 0), c(2, 1)))
#'
#' @export
var_sim_tvp <- function(
    n,
    sigma,
    phi_list,
    T_per_period = NULL,
    mat_tvp = NULL,
    type = c("piecewise", "continuous"),
    burn = 500,
    intercept = 0
) {

  type <- match.arg(type)
  d <- nrow(phi_list[[1]])

  if (type == "piecewise") {
    # Piecewise constant TVP
    n_periods <- length(phi_list)

    if (is.null(T_per_period)) {
      T_per_period <- rep(floor(n / n_periods), n_periods)
      T_per_period[n_periods] <- n - sum(T_per_period[-n_periods])
    }

    stopifnot(length(T_per_period) == n_periods)
    stopifnot(sum(T_per_period) == n)

    # Handle intercept: convert to list of vectors (one per period)
    if (is.list(intercept)) {
      # Already a list - should be one vector per period
      stopifnot(length(intercept) == n_periods)
      intercept_list <- intercept
    } else if (is.numeric(intercept) && length(intercept) == d) {
      # Vector - replicate for each period
      intercept_list <- replicate(n_periods, intercept, simplify = FALSE)
    } else {
      # Scalar - replicate for each period
      intercept_list <- replicate(n_periods, rep(intercept, d), simplify = FALSE)
    }

    # Expand phi and intercept to one per timepoint
    phi_per_timepoint <- unlist(lapply(seq_along(phi_list), function(p) {
      replicate(T_per_period[p], phi_list[[p]], simplify = FALSE)
    }), recursive = FALSE)

    intercept_per_timepoint <- unlist(lapply(seq_along(intercept_list), function(p) {
      replicate(T_per_period[p], intercept_list[[p]], simplify = FALSE)
    }), recursive = FALSE)

    # Simulate using var_sim_growth
    data <- var_sim_growth(
      n = n,
      phi = phi_per_timepoint,
      sigma = sigma,
      burn = burn,
      intercept = intercept_per_timepoint
    )
    colnames(data) <- paste0("V", 1:d)

    # Return structured output
    list(
      data = list(data),
      mat_ind_final = list(phi_list),
      T_per_period = T_per_period,
      breaks = cumsum(T_per_period[-n_periods]),
      intercept = intercept_list
    )

  } else {
    # Continuous (smooth) TVP
    stopifnot(!is.null(mat_tvp))
    stopifnot(length(phi_list) == 1)

    # For continuous, intercept is constant over time
    if (is.numeric(intercept) && length(intercept) == 1) {
      intercept <- rep(intercept, d)
    }

    phi_base <- phi_list[[1]]

    # Generate smoothly varying phi for each timepoint
    phi_per_timepoint <- replicate(n + burn, phi_base, simplify = FALSE)

    for (j in 1:d) {
      for (m in 1:d) {
        if (mat_tvp[j, m]) {
          tvp_series <- cumsum(sample(c(-1, 1), (n + burn), TRUE))
          tvp_series <- smooth(scales::rescale(tvp_series, to = c(-1, 1)))

          phi_per_timepoint <- lapply(seq_along(phi_per_timepoint), function(v) {
            phi_per_timepoint[[v]][j, m] <- tvp_series[v]
            phi_per_timepoint[[v]]
          })
        }
      }
    }

    # Simulate using the phi sequence
    k <- d
    p <- 1
    inno <- MASS::mvrnorm(n = n + burn, rep(0, k), sigma)
    init <- MASS::mvrnorm(n = p, rep(0, k), sigma)
    init <- matrix(init, nrow = p)
    id <- seq(from = p, to = 1, by = -1)
    Y <- matrix(0, (n + burn), k)

    for (r in 1:(n + burn)) {
      Y[r, ] <- intercept + phi_per_timepoint[[r]] %*% as.vector(t(init[id, , drop = FALSE])) + inno[r, ]
      init <- rbind(init[-1, , drop = FALSE], Y[r, ])
    }

    data <- Y[-(1:burn), ]
    colnames(data) <- paste0("V", 1:d)
    phi_final <- phi_per_timepoint[-(1:burn)]

    # Return structured output
    list(
      data = list(data),
      mat_ind_final = list(phi_final),
      intercept = intercept
    )
  }
}

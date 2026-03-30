#' Simulate multivar TVP data (piecewise or continuous)
#'
#' @param k Integer. Number of subjects.
#' @param d Integer. Number of variables per subject.
#' @param n Integer. Total time series length per subject.
#' @param sigma Matrix. Innovation covariance matrix (d x d).
#' @param phi_common_list List of d x d common transition matrices. For piecewise: one per period.
#'   For continuous: one base matrix.
#' @param phi_unique_list List of lists. phi_unique_list[[subject]][[period]] gives the unique
#'   effects matrix for that subject/period. For continuous: phi_unique_list[[subject]] is one matrix.
#' @param T_per_period Integer vector. Observations per period (piecewise only).
#'   If NULL, periods are equal length.
#' @param mat_tvp Logical matrix (d x d). Which elements are time-varying (continuous only).
#' @param type Character. Either "piecewise" or "continuous".
#' @param burn Integer. Burn-in period to discard. Default 500.
#' @param intercept Intercept term. Can be:
#'   - NULL (zeros for all)
#'   - List of vectors: intercept[[subject]] is length d (constant over periods)
#'   - List of lists: intercept[[subject]][[period]] is length d (period-specific)
#'
#' @return List with components:
#'   \item{data}{List of T x d data matrices, one per subject}
#'   \item{mat_com}{Common transition matrix (piecewise: list per period; continuous: list per timepoint)}
#'   \item{mat_ind_unique}{List of unique effects (subject x period)}
#'   \item{mat_ind_final}{List of total effects (subject x period)}
#'   \item{T_per_period}{Observations per period}
#'   \item{breaks}{Break points for period changes}
#'
#' @examples
#' # Piecewise TVP with 2 subjects, 2 periods
#' d <- 3
#' phi_com1 <- matrix(0, d, d); phi_com1[1,2] <- 0.3
#' phi_com2 <- matrix(0, d, d); phi_com2[1,2] <- 0.3; phi_com2[2,3] <- 0.2
#'
#' phi_uniq <- list(
#'   list(matrix(0, d, d), matrix(0, d, d)),  # subject 1: no unique effects
#'   list(matrix(0, d, d), matrix(0, d, d))   # subject 2: no unique effects
#' )
#' phi_uniq[[2]][[1]][3,1] <- 0.2  # subject 2, period 1 has unique edge
#'
#' sim <- multivar_sim_tvp(
#'   k = 2, d = d, n = 200, sigma = diag(d),
#'   phi_common_list = list(phi_com1, phi_com2),
#'   phi_unique_list = phi_uniq,
#'   T_per_period = c(100, 100),
#'   type = "piecewise"
#' )
#'
#' @export
multivar_sim_tvp <- function(
    k,
    d,
    n,
    sigma,
    phi_common_list,
    phi_unique_list = NULL,
    T_per_period = NULL,
    mat_tvp = NULL,
    type = c("piecewise", "continuous"),
    burn = 500,
    intercept = NULL
) {

  type <- match.arg(type)

  if (type == "piecewise") {
    # Piecewise constant TVP
    n_periods <- length(phi_common_list)

    if (is.null(T_per_period)) {
      T_per_period <- rep(floor(n / n_periods), n_periods)
      T_per_period[n_periods] <- n - sum(T_per_period[-n_periods])
    }

    stopifnot(length(T_per_period) == n_periods)
    stopifnot(sum(T_per_period) == n)

    # Handle intercept: convert to list of lists (subject x period)
    if (is.null(intercept)) {
      # Default: zeros for all
      intercept_list <- lapply(1:k, function(i) {
        replicate(n_periods, rep(0, d), simplify = FALSE)
      })
    } else if (is.list(intercept) && !is.list(intercept[[1]])) {
      # List of vectors: replicate each subject's intercept for all periods
      intercept_list <- lapply(intercept, function(int_i) {
        replicate(n_periods, int_i, simplify = FALSE)
      })
    } else {
      # Already list of lists: intercept[[subject]][[period]]
      intercept_list <- intercept
    }

    # Default unique effects: all zeros
    if (is.null(phi_unique_list)) {
      phi_unique_list <- lapply(1:k, function(i) {
        lapply(1:n_periods, function(p) matrix(0, d, d))
      })
    }

    # Build total effects: common + unique for each subject/period
    phi_total_list <- lapply(1:k, function(i) {
      lapply(1:n_periods, function(p) {
        phi_common_list[[p]] + phi_unique_list[[i]][[p]]
      })
    })

    # Simulate data for each subject
    data_list <- lapply(1:k, function(i) {
      # Expand phi to one per timepoint for this subject
      phi_per_timepoint <- unlist(lapply(1:n_periods, function(p) {
        replicate(T_per_period[p], phi_total_list[[i]][[p]], simplify = FALSE)
      }), recursive = FALSE)

      # Expand intercept to one per timepoint for this subject
      intercept_per_timepoint <- unlist(lapply(1:n_periods, function(p) {
        replicate(T_per_period[p], intercept_list[[i]][[p]], simplify = FALSE)
      }), recursive = FALSE)

      # Simulate using var_sim_growth
      dat <- var_sim_growth(
        n = n,
        phi = phi_per_timepoint,
        sigma = sigma,
        burn = burn,
        intercept = intercept_per_timepoint
      )
      colnames(dat) <- paste0("V", 1:d)
      dat
    })

    intercept <- intercept_list

    breaks <- cumsum(T_per_period[-n_periods])

    # Return structured output
    list(
      data = data_list,
      mat_com = phi_common_list,
      mat_ind_unique = phi_unique_list,
      mat_ind_final = phi_total_list,
      T_per_period = T_per_period,
      breaks = breaks,
      intercept = intercept
    )

  } else {
    # Continuous (smooth) TVP
    stopifnot(!is.null(mat_tvp))
    stopifnot(length(phi_common_list) == 1)

    phi_base <- phi_common_list[[1]]

    # Handle intercept: for continuous, constant over time
    if (is.null(intercept)) {
      intercept <- replicate(k, rep(0, d), simplify = FALSE)
    }

    # Default unique effects: all zeros (one matrix per subject)
    if (is.null(phi_unique_list)) {
      phi_unique_list <- lapply(1:k, function(i) matrix(0, d, d))
    }

    # Generate smoothly varying common component
    phi_common_per_timepoint <- replicate(n + burn, phi_base, simplify = FALSE)

    for (j in 1:d) {
      for (m in 1:d) {
        if (mat_tvp[j, m]) {
          tvp_series <- cumsum(sample(c(-1, 1), (n + burn), TRUE))
          tvp_series <- smooth(scales::rescale(tvp_series, to = c(-1, 1)))

          phi_common_per_timepoint <- lapply(seq_along(phi_common_per_timepoint), function(v) {
            phi_common_per_timepoint[[v]][j, m] <- tvp_series[v]
            phi_common_per_timepoint[[v]]
          })
        }
      }
    }

    # Simulate data for each subject
    data_list <- lapply(1:k, function(i) {
      # Total = common + unique (unique is constant over time)
      phi_per_timepoint <- lapply(phi_common_per_timepoint, function(phi_com) {
        phi_com + phi_unique_list[[i]]
      })

      # Simulate
      p <- 1
      inno <- MASS::mvrnorm(n = n + burn, rep(0, d), sigma)
      init <- MASS::mvrnorm(n = p, rep(0, d), sigma)
      init <- matrix(init, nrow = p)
      id <- seq(from = p, to = 1, by = -1)
      Y <- matrix(0, (n + burn), d)

      int_i <- intercept[[i]]

      for (r in 1:(n + burn)) {
        Y[r, ] <- int_i + phi_per_timepoint[[r]] %*% as.vector(t(init[id, , drop = FALSE])) + inno[r, ]
        init <- rbind(init[-1, , drop = FALSE], Y[r, ])
      }

      dat <- Y[-(1:burn), ]
      colnames(dat) <- paste0("V", 1:d)
      dat
    })

    # Extract post-burn-in phi sequences
    phi_common_final <- phi_common_per_timepoint[-(1:burn)]

    # Total effects for each subject (list of lists)
    phi_total_list <- lapply(1:k, function(i) {
      lapply(phi_common_final, function(phi_com) {
        phi_com + phi_unique_list[[i]]
      })
    })

    # Return structured output
    list(
      data = data_list,
      mat_com = phi_common_final,
      mat_ind_unique = phi_unique_list,
      mat_ind_final = phi_total_list,
      intercept = intercept
    )
  }
}

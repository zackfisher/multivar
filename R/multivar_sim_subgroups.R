#' Simulate multivar data with explicit subgroup structure
#'
#' Generates VAR data with a hierarchical structure: common effects (shared by all),
#' subgroup effects (shared within groups), and unique effects (individual-specific).
#' Importantly, these three components occupy NON-OVERLAPPING positions in the
#' coefficient matrices, ensuring clean identifiability.
#'
#' @param k Total number of subjects (sum across all subgroups)
#' @param d Number of variables
#' @param n Number of timepoints
#' @param subgroup Integer vector of length k indicating subgroup membership.
#'   For example, c(1,1,1,2,2,2) indicates 3 subjects in group 1 and 3 in group 2.
#' @param p_com Proportion of positions allocated to common effects.
#'   Sampled from diagonal positions. Default: d/d^2 (i.e., diagonal only).
#' @param p_sub Proportion of positions allocated to subgroup effects (per group).
#'   Sampled from off-diagonal positions excluding common. Default: 0.05.
#' @param p_ind Proportion of positions allocated to unique effects (per individual).
#'   Sampled from remaining positions. Default: 0.05.
#' @param sigma Innovation covariance matrix. If NULL, uses identity matrix.
#' @param lb Lower bound for coefficient values. Default: 0.
#' @param ub Upper bound for coefficient values. Default: 1.
#' @param intercept Optional list of intercept vectors (one per subject).
#'   If NULL, no intercepts are used. Default: NULL.
#' @param max_attempts Maximum number of attempts to generate stationary dynamics.
#'   Default: 100.
#' @param max_eig_threshold Maximum absolute eigenvalue allowed for stationarity.
#'   Default: 0.99.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{data}{List of k data matrices (each n x d)}
#'     \item{subgroup}{Vector of subgroup memberships}
#'     \item{mat_com}{Common effect matrix (d x d)}
#'     \item{mat_sub_unique}{List of k subgroup effect matrices (replicated within groups)}
#'     \item{mat_ind_unique}{List of k unique effect matrices}
#'     \item{mat_ind_final}{List of k total effect matrices}
#'     \item{intercept}{List of intercept vectors (if provided)}
#'   }
#'
#' @details
#' The function generates a three-level hierarchical structure where:
#' \describe{
#'   \item{Common effects}{Shared by all subjects, typically on diagonal positions}
#'   \item{Subgroup effects}{Shared within subgroups, on distinct off-diagonal positions}
#'   \item{Unique effects}{Individual-specific, on remaining positions}
#' }
#'
#' The key feature is that these three components occupy NON-OVERLAPPING positions
#' in the coefficient matrices, which ensures identifiability of the decomposition:
#' \deqn{total[i] = common + subgroup[g] + unique[i]}
#'
#' The function repeatedly generates matrices until all subjects have stationary
#' dynamics (maximum absolute eigenvalue < max_eig_threshold).
#'
#' @examples
#' \dontrun{
#' # Simulate 6 subjects in 2 subgroups
#' sim <- multivar_sim_subgroups(
#'   k = 6,
#'   d = 4,
#'   n = 100,
#'   subgroup = c(1, 1, 1, 2, 2, 2),
#'   p_com = 0.2,
#'   p_sub = 0.1,
#'   p_ind = 0.05
#' )
#'
#' # Fit model
#' object <- constructModel(sim$data, subgroup_membership = sim$subgroup)
#' fit <- cv.multivar(object)
#'
#' # Evaluate including subgroup effects
#' perf <- eval_multivar_performance(sim, fit)
#' print(perf)
#' }
#'
#' @export
multivar_sim_subgroups <- function(k,
                                   d,
                                   n,
                                   subgroup,
                                   p_com = d / d^2,
                                   p_sub = 0.05,
                                   p_ind = 0.05,
                                   sigma = NULL,
                                   lb = 0,
                                   ub = 1,
                                   intercept = NULL,
                                   max_attempts = 100,
                                   max_eig_threshold = 0.99) {

  # Validate inputs
  if (length(subgroup) != k) {
    stop("Length of subgroup must equal k")
  }

  if (any(subgroup < 1)) {
    stop("Subgroup memberships must be positive integers")
  }

  s <- length(unique(subgroup))  # Number of subgroups

  if (is.null(sigma)) {
    sigma <- diag(d)
  }

  if (nrow(sigma) != d || ncol(sigma) != d) {
    stop("sigma must be a d x d matrix")
  }

  # Attempt to generate stationary system
  nonstationary <- TRUE
  attempts <- 0

  while (nonstationary && attempts < max_attempts) {
    attempts <- attempts + 1

    # Set up true matrices
    true_com <- matrix(0, d, d, byrow = TRUE)
    true_sub <- array(0, dim = c(d, d, k))
    true_ind <- array(0, dim = c(d, d, k))

    # Set up position indices
    diag_pos <- 1 + 0:(d - 1) * (d + 1)

    # Common positions: sample from diagonal
    n_com <- round(p_com * d^2)
    if (n_com > length(diag_pos)) {
      n_com <- length(diag_pos)
    }
    com_pos <- if (n_com > 0) sample(x = diag_pos, size = n_com) else integer(0)

    # Subgroup positions: sample from off-diagonal, excluding common
    n_sub <- round(p_sub * s * d^2)
    available_for_sub <- setdiff(1:d^2, c(diag_pos, com_pos))

    if (n_sub > length(available_for_sub)) {
      n_sub <- length(available_for_sub)
    }

    sub_pos <- if (n_sub > 0) {
      sample(x = available_for_sub, size = n_sub)
    } else {
      integer(0)
    }

    # Split subgroup positions among s subgroups
    if (length(sub_pos) > 0) {
      sub_pos_by_group <- split(sub_pos, factor(sort(rank(sub_pos) %% s)))
    } else {
      sub_pos_by_group <- rep(list(integer(0)), s)
    }

    # Replicate subgroup positions to all subjects based on membership
    sub_pos_by_subject <- lapply(1:k, function(i) {
      sub_pos_by_group[[subgroup[i]]]
    })

    # Unique positions: sample from remaining positions (per subject)
    n_ind <- round(p_ind * d^2)
    ind_pos <- lapply(seq_len(k), function(i) {
      available_for_ind <- setdiff(1:d^2, c(diag_pos, com_pos, unlist(sub_pos_by_subject)))
      if (n_ind > length(available_for_ind)) {
        n_ind_i <- length(available_for_ind)
      } else {
        n_ind_i <- n_ind
      }

      if (n_ind_i > 0) {
        sample(x = available_for_ind, size = n_ind_i, replace = FALSE)
      } else {
        integer(0)
      }
    })

    # Fill matrices
    if (length(com_pos) > 0) {
      true_com[com_pos] <- runif(length(com_pos), lb, ub)
    }

    # Generate subgroup matrices (one per subgroup)
    true_sub_by_group <- lapply(seq_len(s), function(g) {
      mat <- matrix(0, d, d)
      if (length(sub_pos_by_group[[g]]) > 0) {
        mat[sub_pos_by_group[[g]]] <- runif(length(sub_pos_by_group[[g]]), lb, ub)
      }
      mat
    })

    # Replicate subgroup matrices to subjects based on membership
    true_sub <- lapply(seq_len(k), function(j) {
      true_sub_by_group[[subgroup[j]]]
    })

    true_ind <- lapply(seq_len(k), function(j) {
      mat <- matrix(0, d, d)
      if (length(ind_pos[[j]]) > 0) {
        mat[ind_pos[[j]]] <- runif(length(ind_pos[[j]]), lb, ub)
      }
      mat
    })

    # Compute total matrices
    mats <- lapply(seq_len(k), function(j) {
      true_com + true_sub[[j]] + true_ind[[j]]
    })

    # Check stationarity (all eigenvalues < threshold)
    nonstationary <- any(unlist(lapply(mats, function(x) {
      max_eig <- max(abs(eigen(x, only.values = TRUE)$values))
      max_eig > max_eig_threshold
    })))
  }

  if (nonstationary) {
    warning(sprintf(
      "Failed to generate stationary system after %d attempts. Returning last attempt (may be nonstationary).",
      max_attempts
    ))
  }

  # Simulate time series data using multivar_sim with mat_total
  sim <- multivar_sim(
    k = k,
    d = d,
    n = n,
    mat_total = mats,
    sigma = sigma,
    intercept = intercept
  )

  # Add ground truth components to output
  sim$subgroup <- subgroup
  sim$mat_com <- true_com
  sim$mat_sub_unique <- true_sub
  sim$mat_ind_unique <- true_ind
  sim$mat_ind_final <- mats

  sim
}

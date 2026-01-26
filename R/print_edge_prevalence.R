#' Print edge prevalence across subjects
#'
#' For models with multiple subjects (K>1), displays a matrix showing
#' the proportion/percentage of subjects that have each edge (non-zero entry).
#'
#' @param fit A fitted multivar object (output from cv.multivar())
#' @param type Character; "proportion" (default) or "count" to show counts instead
#' @param digits Number of decimal places (default: 2)
#' @param zero_char Character to display for edges present in 0 subjects (default: ".")
#' @param threshold Threshold for considering an edge as present (default: 1e-10)
#'
#' @return Invisibly returns the prevalence matrix
#' @export
#'
#' @examples
#' # fit <- cv.multivar(object)  # K=2 or more subjects
#' # print_edge_prevalence(fit)
print_edge_prevalence <- function(fit,
                                  type = c("proportion", "count"),
                                  digits = 2,
                                  zero_char = ".",
                                  threshold = 1e-10) {

  type <- match.arg(type)

  # Extract object
  obj <- fit$obj
  if (is.null(obj)) {
    stop("Fit object does not contain the original model specification (missing $obj)")
  }

  k <- obj@k
  tvp <- obj@tvp

  # Check if K > 1
  if (k <= 1) {
    message("Edge prevalence is only meaningful for K>1 models (multiple subjects).")
    message("For K=1 models, use print_dynamics(fit) instead.")
    return(invisible(NULL))
  }

  # Get total effects
  total_mats <- fit$mats$total

  # Header
  cat("\n")
  cat(strrep("=", 70), "\n")
  if (type == "proportion") {
    cat("Edge Prevalence (Proportion of Subjects)\n")
  } else {
    cat("Edge Prevalence (Count of Subjects)\n")
  }
  cat(strrep("=", 70), "\n\n")

  if (tvp) {
    # TVP model: show prevalence for each period
    n_periods <- length(total_mats[[1]])

    for (p in seq_len(n_periods)) {
      cat(sprintf("Period %d:\n", p))
      cat(strrep("-", 70), "\n")

      # Extract period p for all subjects
      period_mats <- lapply(total_mats, function(x) x[[p]])

      # Compute prevalence matrix
      prev_mat <- compute_prevalence(period_mats, k, threshold)

      # Format and print
      print_prevalence_matrix(prev_mat, type, digits, zero_char, k)
      cat("\n")
    }
  } else {
    # Non-TVP model
    # Compute prevalence matrix
    prev_mat <- compute_prevalence(total_mats, k, threshold)

    # Format and print
    print_prevalence_matrix(prev_mat, type, digits, zero_char, k)
  }

  cat(strrep("=", 70), "\n\n")

  invisible(prev_mat)
}


# Helper function to compute prevalence matrix
compute_prevalence <- function(mat_list, k, threshold) {
  # Get dimensions from first matrix
  d <- nrow(mat_list[[1]])

  # Count how many subjects have each edge
  count_mat <- matrix(0, nrow = d, ncol = d)

  for (i in seq_len(k)) {
    mat <- mat_list[[i]]
    count_mat <- count_mat + (abs(mat) > threshold)
  }

  # Convert to proportion
  prev_mat <- count_mat / k

  # Preserve row/column names if available
  if (!is.null(rownames(mat_list[[1]]))) {
    rownames(prev_mat) <- rownames(mat_list[[1]])
  } else {
    rownames(prev_mat) <- paste0("V", 1:d, "_t")
  }

  if (!is.null(colnames(mat_list[[1]]))) {
    colnames(prev_mat) <- colnames(mat_list[[1]])
  } else {
    colnames(prev_mat) <- paste0("V", 1:d, "_{t-1}")
  }

  return(list(proportion = prev_mat, count = count_mat))
}


# Helper function to print prevalence matrix
print_prevalence_matrix <- function(prev_list, type, digits, zero_char, k) {
  if (type == "proportion") {
    mat_to_print <- prev_list$proportion
  } else {
    mat_to_print <- prev_list$count
  }

  d <- nrow(mat_to_print)

  # Format values
  mat_formatted <- matrix("", nrow = d, ncol = d)
  for (i in 1:d) {
    for (j in 1:d) {
      val <- mat_to_print[i, j]
      if (val == 0) {
        mat_formatted[i, j] <- zero_char
      } else {
        if (type == "proportion") {
          # Show as proportion (0.00 to 1.00)
          mat_formatted[i, j] <- sprintf(paste0("%.", digits, "f"), val)
        } else {
          # Show as count
          mat_formatted[i, j] <- sprintf("%d/%d", round(val), k)
        }
      }
    }
  }

  # Add row and column names
  rownames(mat_formatted) <- rownames(mat_to_print)
  colnames(mat_formatted) <- colnames(mat_to_print)

  # Print
  print(noquote(mat_formatted))

  # Summary statistics
  cat("\n")
  if (type == "proportion") {
    cat(sprintf("  Universal edges (100%%):        %d\n",
                sum(prev_list$proportion == 1.0)))
    cat(sprintf("  Common edges (>50%%):           %d\n",
                sum(prev_list$proportion > 0.5)))
    cat(sprintf("  Rare edges (<50%%):             %d\n",
                sum(prev_list$proportion > 0 & prev_list$proportion < 0.5)))
    cat(sprintf("  Absent edges (0%%):             %d\n",
                sum(prev_list$proportion == 0)))
  } else {
    cat(sprintf("  Universal edges (all %d):       %d\n",
                k, sum(prev_list$count == k)))
    cat(sprintf("  Variable edges (between 1-%d):  %d\n",
                k-1, sum(prev_list$count > 0 & prev_list$count < k)))
    cat(sprintf("  Absent edges (0):               %d\n",
                sum(prev_list$count == 0)))
  }
}

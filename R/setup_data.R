#' @importFrom Matrix Matrix
#' @export
setup_data <- function (data, standardize, lag, horizon, intercept, tvp = FALSE, breaks = NULL) {
  
   if (is.null(data)){
    stop(paste0(
      "multivar ERROR: neither a data directory nor a data list is specified. ",
      "Please either specify a directory of data or a list of individual data files."
    )
    )
  }
  
  #-------------------------------------------------------------#
  # If the data is already in list form.
  #-------------------------------------------------------------#
  if (!is.list(data)){
      stop(paste0(
        "multivar ERROR: data must be supplied as a list of matrices."
      ))
  } else {
      ts_list  <- list()
  
      # if the user-supplied list does not have names, add names
      if(is.null(names(data))){
  
        names(data) <- paste0("dataset", 1:length(data))
  
      }
  
      ts_list  <- data
  }

  
  #-------------------------------------------------------------#
  # Ensure all datafiles share the same column order.
  #
  #   # Now we also make all datasets have the same column order.
  #-------------------------------------------------------------# 
  varnames       <- colnames(ts_list[[1]])
  n_orig_vars    <- ncol(ts_list[[1]])
  
  if(is.null(varnames)){
    varnames <- c(paste0("V", seq(1,n_orig_vars)))
    ts_list <- lapply(ts_list, function(x) { colnames(x) <- varnames;  x })
  } 
  
  
  #-------------------------------------------------------------#
  # Final data manipulation.
  #-------------------------------------------------------------#
  if (ncol(ts_list[[1]]) == 1) {
    stop(paste0("multivar ERROR: only one column of data read in. ",
                "Check if sep argument properly specified."))
  }
  
  

  # For non-TVP: global standardization on raw time series
  # For TVP: defer standardization until after A/b split (per-period)
  # Note: Always center data for good LASSO estimation. Intercepts are recovered post-hoc.
  if(standardize & !tvp){
    ts_list <- lapply(ts_list, function(df) {scale(df)})
  }


  #-------------------------------------------------------------#
  # Final data checks
  #-------------------------------------------------------------#
  
  n_subjects   <- length(ts_list)
  cols         <- numeric()
  missingData  <- numeric()
  constantCols <- logical()
  numericCols  <- logical()
  
  # check for obvious errors in data
  for (k in 1:length(ts_list)){
    data.file <- ts_list[[k]]
    cols[k]   <- ncol(data.file)
    missingData[k]  <- any(is.na(data.file))
    constantCols[k] <- any(apply(data.file, 2, sd, na.rm = TRUE) == 0)
    numericCols[k]  <- any(apply(data.file, 2, is.numeric) == FALSE)
  }
  
  # Helper function for splitting by breaks

  splitAt <- function(x, pos) {
    unname(split(x, cumsum(seq_along(x) %in% pos)))
  }

  # will need to be updated for lags > 1
  ts_list <- lapply(seq_along(ts_list), function(k){
    df <- ts_list[[k]]
    if(horizon > 0){
      H  <- df[((nrow(df)-horizon+1):nrow(df)),, drop = FALSE]
      df <- df[-((nrow(df)-horizon+1):nrow(df)),, drop = FALSE]
    } else {
      H  <- NA
    }
    A  <- Matrix(df[1:(nrow(df)-1), ], sparse = FALSE)
    b  <- Matrix(df[2:(nrow(df)  ), ], sparse = FALSE)
    colnames(A) <- colnames(b) <- colnames(df)

    # Store original means for intercept recovery (before any standardization)
    # These are needed for: c = mean(b) - Phi * mean(A)
    mean_A_orig <- colMeans(A)
    mean_b_orig <- colMeans(b)

    # Initialize storage for standardization parameters and period means
    sd_A <- sd_b <- NULL
    mean_A_periods <- mean_b_periods <- NULL
    sd_A_periods <- sd_b_periods <- NULL

    # Per-period standardization for TVP
    # Key: compute standardization parameters from A only, apply to both A and b
    # This preserves the VAR relationship b = Phi*A + epsilon
    if(standardize & tvp & !is.null(breaks)){
      ntk <- nrow(A)
      period_indices <- splitAt(seq_len(ntk), breaks[[k]])

      A_standardized <- A
      b_standardized <- b

      # Storage for per-period means and SDs (for intercept recovery)
      mean_A_periods <- vector("list", length(period_indices))
      mean_b_periods <- vector("list", length(period_indices))
      sd_A_periods <- vector("list", length(period_indices))
      sd_b_periods <- vector("list", length(period_indices))

      for(p in seq_along(period_indices)){
        idx <- period_indices[[p]]
        A_period <- A[idx, , drop = FALSE]
        b_period <- b[idx, , drop = FALSE]

        # Compute parameters from A only
        col_means_A <- colMeans(A_period)
        col_means_b <- colMeans(b_period)
        col_sds <- apply(A_period, 2, sd)
        col_sds[col_sds == 0] <- 1  # Avoid division by zero

        # Store period means and SDs for intercept recovery
        mean_A_periods[[p]] <- as.vector(col_means_A)
        mean_b_periods[[p]] <- as.vector(col_means_b)
        sd_A_periods[[p]] <- as.vector(col_sds)
        sd_b_periods[[p]] <- as.vector(apply(b_period, 2, sd))
        sd_b_periods[[p]][sd_b_periods[[p]] == 0] <- 1

        # Always center and scale for good LASSO estimation
        # Intercepts are recovered post-hoc using stored means
        A_standardized[idx, ] <- scale(A_period, center = col_means_A, scale = col_sds)
        b_standardized[idx, ] <- scale(b_period, center = col_means_A, scale = col_sds)
      }

      A <- A_standardized
      b <- b_standardized

    } else if (tvp & !is.null(breaks)) {
      # TVP without standardization: still need period means for intercept recovery
      ntk <- nrow(A)
      period_indices <- splitAt(seq_len(ntk), breaks[[k]])

      mean_A_periods <- vector("list", length(period_indices))
      mean_b_periods <- vector("list", length(period_indices))

      for(p in seq_along(period_indices)){
        idx <- period_indices[[p]]
        A_period <- A[idx, , drop = FALSE]
        b_period <- b[idx, , drop = FALSE]

        mean_A_periods[[p]] <- as.vector(colMeans(A_period))
        mean_b_periods[[p]] <- as.vector(colMeans(b_period))
      }
      # No sd_*_periods needed when not standardizing

    } else if (standardize & !tvp) {
      # Non-TVP standardization: store global SDs for intercept recovery
      sd_A <- apply(A, 2, sd)
      sd_b <- apply(b, 2, sd)
      sd_A[sd_A == 0] <- 1
      sd_b[sd_b == 0] <- 1
    }

    # Return data with means for intercept recovery
    list(
      b = b,
      A = A,
      H = H,
      mean_A = as.vector(mean_A_orig),
      mean_b = as.vector(mean_b_orig),
      sd_A = sd_A,
      sd_b = sd_b,
      mean_A_periods = mean_A_periods,
      mean_b_periods = mean_b_periods,
      sd_A_periods = sd_A_periods,
      sd_b_periods = sd_b_periods
    )
  })

  # Note: Intercept column is NOT added to design matrix.
  # Intercepts are recovered post-hoc using: c = mean(b) - Phi * mean(A)
  # See recover_intercepts.R

  return(ts_list)
  
}
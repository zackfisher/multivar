#' Summarize multivar performance against simulation truth (robust to missing parts)
#'
#' Compare estimated transition matrices to simulation truth, but only for
#' components that actually exist in both `sim_obj` and `fit_obj`. For example,
#' if `fit_obj$mats` only contains `total`, only "total" rows are returned.
#'
#' @param sim_obj  Output of multivar_sim()
#' @param fit_obj  Output of cv.multivar() / cv.fused() (or similar list with $mats)
#' @param eps      Small positive constant to stabilize relative errors
#' @param reduced.output If TRUE, keep a reduced column set
#' @param averages.only  If TRUE, return only summary rows (unique_mean, total_mean)
#'                       and the single 'common' row (if present)
#' @param label Optional label for the fitted model; by default, uses the
#'                   deparsed name of `fit_obj`
#' @param intercept Logical; if TRUE, evaluates intercepts separately from dynamics
#' @export
eval_multivar_performance <- function(sim_obj,
                                      fit_obj,
                                      eps = 1e-8,
                                      reduced.output = TRUE,
                                      averages.only = TRUE,
                                      label = NULL,
                                      intercept = FALSE) {
  ## derive a label for the fitted object
  model_label <- if (is.null(label)) {
    # try to capture the symbol used at call site
    lb <- try(deparse(substitute(fit_obj)), silent = TRUE)
    if (inherits(lb, "try-error")) "fit_obj" else lb
  } else {
    as.character(label)
  }
  
  ## ---- helpers ----
  as_list <- function(x) if (is.null(x)) list() else if (is.list(x)) x else list(x)
  is_mat  <- function(x) is.matrix(x) && is.numeric(x)

  ## ---- intercept helpers ----
  has_intercept <- function(mat) {
    is_mat(mat) && !is.null(colnames(mat)) && colnames(mat)[1] == "Intercept"
  }

  extract_intercept <- function(mat) {
    if (has_intercept(mat)) mat[, 1, drop = TRUE] else NULL
  }

  extract_dynamics <- function(mat) {
    if (has_intercept(mat)) mat[, -1, drop = FALSE] else mat
  }

  compute_metrics_for_vector <- function(true_vec, est_vec, effect_label,
                                         subject_id = NA_integer_,
                                         eps = 1e-8, msfe_min_global = NA_real_) {
    # For intercepts (vectors), compute the same metrics as for matrices
    # Treat vector elements as edges for detection metrics

    # Binary edge indicators
    true_nz <- (abs(true_vec) > 0)
    est_nz  <- (abs(est_vec)  > 0)

    TP <- sum(true_nz & est_nz)
    TN <- sum(!true_nz & !est_nz)
    FP <- sum(!true_nz &  est_nz)
    FN <- sum( true_nz & !est_nz)

    sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
    prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    f1   <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0) 2 * prec * sens / (prec + sens) else NA_real_

    denom_mcc <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    mcc <- if (denom_mcc > 0) ((TP * TN) - (FP * FN)) / denom_mcc else NA_real_

    diff_vec <- est_vec - true_vec

    # Overall bias metrics
    abs_bias <- mean(abs(diff_vec))
    rel_bias <- mean(abs(diff_vec) / (abs(true_vec) + eps))
    rmse     <- sqrt(mean(diff_vec^2))
    mae      <- mean(abs(diff_vec))

    # Variability over all entries
    sd_error <- sqrt(mean((diff_vec - mean(diff_vec))^2))

    # Bias & variability restricted to truly nonzero entries
    if (any(true_nz)) {
      diff_nz       <- diff_vec[true_nz]
      true_nz_vals  <- true_vec[true_nz]
      abs_bias_nz   <- mean(abs(diff_nz))
      rel_bias_nz   <- mean(abs(diff_nz) / (abs(true_nz_vals) + eps))
      sd_error_nz   <- sqrt(mean((diff_nz - mean(diff_nz))^2))
    } else {
      abs_bias_nz <- NA_real_
      rel_bias_nz <- NA_real_
      sd_error_nz <- NA_real_
    }

    data.frame(
      model        = model_label,
      effect       = effect_label,
      subject      = subject_id,
      TP           = TP,
      FP           = FP,
      TN           = TN,
      FN           = FN,
      sensitivity  = sens,
      specificity  = spec,
      precision    = prec,
      F1           = f1,
      MCC          = mcc,
      abs_bias     = abs_bias,
      rel_bias     = rel_bias,
      abs_bias_nz  = abs_bias_nz,
      rel_bias_nz  = rel_bias_nz,
      sd_error     = sd_error,
      sd_error_nz  = sd_error_nz,
      RMSE         = rmse,
      MAE          = mae,
      MSFE_min     = msfe_min_global,
      stringsAsFactors = FALSE
    )
  }
  
  compute_metrics_for_pair <- function(true_mat, est_mat, effect_label,
                                       subject_id = NA_integer_,
                                       eps = 1e-8, msfe_min_global = NA_real_) {
    true_vec <- as.numeric(true_mat)
    est_vec  <- as.numeric(est_mat)
    
    # binary edge indicators
    true_nz <- (abs(true_vec) > 0)
    est_nz  <- (abs(est_vec)  > 0)
    
    TP <- sum(true_nz & est_nz)
    TN <- sum(!true_nz & !est_nz)
    FP <- sum(!true_nz &  est_nz)
    FN <- sum( true_nz & !est_nz)
    
    sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
    prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    f1   <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0) 2 * prec * sens / (prec + sens) else NA_real_
    
    denom_mcc <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    mcc <- if (denom_mcc > 0) ((TP * TN) - (FP * FN)) / denom_mcc else NA_real_
    
    diff_vec <- est_vec - true_vec
    
    ## overall bias metrics (unchanged)
    abs_bias <- mean(abs(diff_vec))
    rel_bias <- mean(abs(diff_vec) / (abs(true_vec) + eps))
    rmse     <- sqrt(mean(diff_vec^2))
    mae      <- mean(abs(diff_vec))
    
    ## NEW: variability over all entries
    sd_error <- sqrt(mean((diff_vec - mean(diff_vec))^2))
    
    ## NEW: bias & variability restricted to truly nonzero entries
    if (any(true_nz)) {
      diff_nz       <- diff_vec[true_nz]
      true_nz_vals  <- true_vec[true_nz]
      abs_bias_nz   <- mean(abs(diff_nz))
      rel_bias_nz   <- mean(abs(diff_nz) / (abs(true_nz_vals) + eps))
      sd_error_nz   <- sqrt(mean((diff_nz - mean(diff_nz))^2))
    } else {
      abs_bias_nz <- NA_real_
      rel_bias_nz <- NA_real_
      sd_error_nz <- NA_real_
    }
    
    data.frame(
      model        = model_label,
      effect       = effect_label,
      subject      = subject_id,
      TP           = TP, FP = FP, TN = TN, FN = FN,
      sensitivity  = sens,
      specificity  = spec,
      precision    = prec,
      F1           = f1,
      MCC          = mcc,
      abs_bias     = abs_bias,
      rel_bias     = rel_bias,
      abs_bias_nz  = abs_bias_nz,
      rel_bias_nz  = rel_bias_nz,
      sd_error     = sd_error,
      sd_error_nz  = sd_error_nz,
      RMSE         = rmse,
      MAE          = mae,
      MSFE_min     = msfe_min_global,
      stringsAsFactors = FALSE
    )
  }
  
  add_summary_block <- function(df, which_effect, label) {
    sub <- df[df$effect == which_effect & !is.na(df$subject), , drop = FALSE]
    if (nrow(sub) == 0) return(NULL)

    # Determine which numeric columns are actually present
    possible_num_cols <- c("TP","FP","TN","FN",
                           "sensitivity","specificity","precision","F1","MCC",
                           "abs_bias","rel_bias","abs_bias_nz","rel_bias_nz",
                           "sd_error","sd_error_nz",
                           "RMSE","MAE","MSFE_min")
    num_cols <- intersect(possible_num_cols, colnames(sub))

    means <- colMeans(sub[, num_cols, drop = FALSE], na.rm = TRUE)
    data.frame(model = model_label,
               effect = label, subject = NA_integer_, t(means),
               stringsAsFactors = FALSE, row.names = NULL)
  }
  
  ## ---- pull truth (optional) ----
  truth_common  <- if (!is.null(sim_obj$mat_com) && is_mat(sim_obj$mat_com)) sim_obj$mat_com else NULL
  truth_uniqueL <- Filter(is_mat, as_list(sim_obj$mat_ind_unique))
  truth_subgrpL <- Filter(is_mat, as_list(sim_obj$mat_sub_unique))
  truth_totalL  <- Filter(is_mat, as_list(sim_obj$mat_ind_final))

  ## ---- detect K=1 case ----
  # For K=1, there's no multi-subject decomposition, so common/unique comparisons are meaningless
  # In K=1, fit outputs common=unique=total, but sim may have separate common/unique components
  # We should only evaluate total effects
  k_subjects <- length(Filter(is_mat, as_list(fit_obj$mats$total)))
  is_k1 <- (k_subjects == 1)

  ## ---- pull estimates (optional) ----
  # For K=1, skip common and unique estimates (they're identical to total and meaningless)
  est_common  <- if (!is_k1 && !is.null(fit_obj$mats$common) && is_mat(fit_obj$mats$common)) fit_obj$mats$common else NULL
  est_uniqueL <- if (!is_k1) Filter(is_mat, as_list(fit_obj$mats$unique)) else list()
  est_subgrpL <- if (!is_k1) Filter(is_mat, as_list(fit_obj$mats$subgrp)) else list()

  ## ---- detect TVP structure (before filtering) ----
  # If fit_obj$mats$total contains lists (not matrices), it's a TVP model
  est_totalL_raw <- as_list(fit_obj$mats$total)
  is_tvp <- length(est_totalL_raw) > 0 && is.list(est_totalL_raw[[1]]) && !is_mat(est_totalL_raw[[1]])

  # Filter based on TVP status
  if (is_tvp) {
    # For TVP: keep list-of-lists structure, don't filter by is_mat
    est_totalL <- est_totalL_raw
    # Also update truth to not filter
    truth_totalL <- as_list(sim_obj$mat_ind_final)

    # Extract TVP components from estimates
    est_common_tvp <- as_list(fit_obj$mats$common_tvp)  # List of matrices (one per period)
    est_unique_tvp <- as_list(fit_obj$mats$tvp)  # List of list of matrices (subjects x time)

    # Extract TVP components from truth (if available)
    truth_common_tvp <- as_list(sim_obj$mat_com_tvp)  # List of matrices (one per period)
    truth_unique_tvp <- as_list(sim_obj$mat_ind_tvp)  # List of list of matrices (subjects x time)
  } else {
    # For non-TVP: filter to keep only matrices
    est_totalL <- Filter(is_mat, est_totalL_raw)
    est_common_tvp <- NULL
    est_unique_tvp <- NULL
    truth_common_tvp <- NULL
    truth_unique_tvp <- NULL
  }

  ## ---- handle intercepts if present ----
  truth_intercepts <- list()
  est_intercepts <- list()

  if (intercept) {
    # Extract intercept vectors from truth
    if (!is.null(sim_obj$intercept)) {
      # True intercepts stored directly in sim_obj$intercept
      truth_intercepts$total <- as_list(sim_obj$intercept)

      # Compute true common intercept (mean across subjects)
      if (length(truth_intercepts$total) > 0) {
        truth_intercepts$common <- rowMeans(do.call(cbind, truth_intercepts$total))
        # Compute true unique intercepts (deviations)
        truth_intercepts$unique <- lapply(truth_intercepts$total,
                                         function(x) x - truth_intercepts$common)
      }
    }

    # Extract intercept vectors from estimates
    if (!is.null(est_common) && has_intercept(est_common)) {
      est_intercepts$common <- extract_intercept(est_common)
    }
    if (length(est_uniqueL) > 0 && has_intercept(est_uniqueL[[1]])) {
      est_intercepts$unique <- lapply(est_uniqueL, extract_intercept)
    }
    if (length(est_totalL) > 0 && has_intercept(est_totalL[[1]])) {
      est_intercepts$total <- lapply(est_totalL, extract_intercept)
    }

    # Extract dynamics (without intercept column) for evaluation
    if (!is.null(est_common) && has_intercept(est_common)) {
      est_common <- extract_dynamics(est_common)
    }
    if (length(est_uniqueL) > 0 && has_intercept(est_uniqueL[[1]])) {
      est_uniqueL <- lapply(est_uniqueL, extract_dynamics)
    }
    if (length(est_totalL) > 0 && has_intercept(est_totalL[[1]])) {
      est_totalL <- lapply(est_totalL, extract_dynamics)
    }

    # Truth doesn't have intercept column in matrices (it's separate)
    # So no need to extract dynamics from truth matrices
  }
  
  ## ---- MSFE (optional) ----
  msfe_min <- NA_real_
  if (!is.null(fit_obj$MSFE)) {
    mm <- try(suppressWarnings(colMeans(fit_obj$MSFE)), silent = TRUE)
    if (!inherits(mm, "try-error")) msfe_min <- suppressWarnings(min(mm, na.rm = TRUE))
  }
  
  out_list <- list()
  
  ## ---- common ----
  if (!is.null(truth_common) && !is.null(est_common)) {
    out_list[[length(out_list) + 1]] <-
      compute_metrics_for_pair(truth_common, est_common, "common",
                               subject_id = NA_integer_, eps = eps,
                               msfe_min_global = msfe_min)
  }
  
  ## ---- unique (match lengths safely) ----
  if (length(truth_uniqueL) > 0 && length(est_uniqueL) > 0) {
    k_u <- min(length(truth_uniqueL), length(est_uniqueL))
    for (ii in seq_len(k_u)) {
      out_list[[length(out_list) + 1]] <-
        compute_metrics_for_pair(truth_uniqueL[[ii]], est_uniqueL[[ii]],
                                 "unique", subject_id = ii, eps = eps,
                                 msfe_min_global = msfe_min)
    }
  }

  ## ---- subgrp (match lengths safely) ----
  if (length(truth_subgrpL) > 0 && length(est_subgrpL) > 0) {
    k_sg <- min(length(truth_subgrpL), length(est_subgrpL))
    for (ii in seq_len(k_sg)) {
      out_list[[length(out_list) + 1]] <-
        compute_metrics_for_pair(truth_subgrpL[[ii]], est_subgrpL[[ii]],
                                 "subgrp", subject_id = ii, eps = eps,
                                 msfe_min_global = msfe_min)
    }
  }

  ## ---- total (match lengths safely) ----
  if (length(truth_totalL) > 0 && length(est_totalL) > 0) {
    k_t <- min(length(truth_totalL), length(est_totalL))

    if (is_tvp) {
      # For TVP models: truth_totalL[[i]] and est_totalL[[i]] are lists of matrices (one per period)
      # Compute metrics by aggregating across time points AND period-specific metrics
      for (ii in seq_len(k_t)) {
        truth_tvp_list <- truth_totalL[[ii]]
        est_tvp_list   <- est_totalL[[ii]]

        # Only proceed if both are lists
        if (is.list(truth_tvp_list) && is.list(est_tvp_list)) {
          n_t <- min(length(truth_tvp_list), length(est_tvp_list))

          if (n_t > 0) {
            # Stack all time-varying matrices into 3D arrays for comparison
            truth_mats <- Filter(is_mat, truth_tvp_list[1:n_t])
            est_mats   <- Filter(is_mat, est_tvp_list[1:n_t])

            if (length(truth_mats) > 0 && length(est_mats) > 0) {
              # Flatten across time to get overall metrics
              truth_vec <- as.numeric(do.call(cbind, truth_mats))
              est_vec   <- as.numeric(do.call(cbind, est_mats))

              # Binary edge indicators
              true_nz <- (abs(truth_vec) > 0)
              est_nz  <- (abs(est_vec)  > 0)

              TP <- sum(true_nz & est_nz)
              TN <- sum(!true_nz & !est_nz)
              FP <- sum(!true_nz &  est_nz)
              FN <- sum( true_nz & !est_nz)

              sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
              spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
              prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
              f1   <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0) 2 * prec * sens / (prec + sens) else NA_real_

              denom_mcc <- sqrt(as.numeric(TP + FP) * as.numeric(TP + FN) * as.numeric(TN + FP) * as.numeric(TN + FN))
              mcc <- if (!is.na(denom_mcc) && denom_mcc > 0) ((TP * TN) - (FP * FN)) / denom_mcc else NA_real_

              diff_vec <- est_vec - truth_vec

              abs_bias <- mean(abs(diff_vec))
              rel_bias <- mean(abs(diff_vec) / (abs(truth_vec) + eps))
              rmse     <- sqrt(mean(diff_vec^2))
              mae      <- mean(abs(diff_vec))
              sd_error <- sqrt(mean((diff_vec - mean(diff_vec))^2))

              if (any(true_nz)) {
                diff_nz       <- diff_vec[true_nz]
                true_nz_vals  <- truth_vec[true_nz]
                abs_bias_nz   <- mean(abs(diff_nz))
                rel_bias_nz   <- mean(abs(diff_nz) / (abs(true_nz_vals) + eps))
                sd_error_nz   <- sqrt(mean((diff_nz - mean(diff_nz))^2))
              } else {
                abs_bias_nz <- NA_real_
                rel_bias_nz <- NA_real_
                sd_error_nz <- NA_real_
              }

              # Overall aggregated metrics
              out_list[[length(out_list) + 1]] <- data.frame(
                model        = model_label,
                effect       = "total",
                subject      = ii,
                TP           = TP, FP = FP, TN = TN, FN = FN,
                sensitivity  = sens,
                specificity  = spec,
                precision    = prec,
                F1           = f1,
                MCC          = mcc,
                abs_bias     = abs_bias,
                rel_bias     = rel_bias,
                abs_bias_nz  = abs_bias_nz,
                rel_bias_nz  = rel_bias_nz,
                sd_error     = sd_error,
                sd_error_nz  = sd_error_nz,
                RMSE         = rmse,
                MAE          = mae,
                MSFE_min     = msfe_min,
                stringsAsFactors = FALSE
              )

              # Period-specific metrics
              for (period_idx in seq_len(n_t)) {
                truth_period <- truth_mats[[period_idx]]
                est_period   <- est_mats[[period_idx]]

                if (is_mat(truth_period) && is_mat(est_period)) {
                  # Compute metrics for this period
                  out_list[[length(out_list) + 1]] <-
                    compute_metrics_for_pair(
                      truth_period, est_period,
                      paste0("total_period", period_idx),
                      subject_id = ii, eps = eps,
                      msfe_min_global = msfe_min
                    )
                }
              }
            }
          }
        }
      }
    } else {
      # Non-TVP case: matrices directly
      for (ii in seq_len(k_t)) {
        out_list[[length(out_list) + 1]] <-
          compute_metrics_for_pair(truth_totalL[[ii]], est_totalL[[ii]],
                                   "total", subject_id = ii, eps = eps,
                                   msfe_min_global = msfe_min)
      }
    }
  }

  ## ---- common_tvp (for TVP models) ----
  if (is_tvp && length(truth_common_tvp) > 0 && length(est_common_tvp) > 0) {
    # Common TVP: list of matrices (one per period)
    # Flatten across periods and evaluate
    n_periods <- min(length(truth_common_tvp), length(est_common_tvp))

    truth_mats_tvp <- Filter(is_mat, truth_common_tvp[1:n_periods])
    est_mats_tvp   <- Filter(is_mat, est_common_tvp[1:n_periods])

    if (length(truth_mats_tvp) > 0 && length(est_mats_tvp) > 0) {
      # Flatten across periods
      truth_vec <- as.numeric(do.call(cbind, truth_mats_tvp))
      est_vec   <- as.numeric(do.call(cbind, est_mats_tvp))

      # Binary edge indicators
      true_nz <- (abs(truth_vec) > 0)
      est_nz  <- (abs(est_vec)  > 0)

      TP <- sum(true_nz & est_nz)
      TN <- sum(!true_nz & !est_nz)
      FP <- sum(!true_nz &  est_nz)
      FN <- sum( true_nz & !est_nz)

      sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
      spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
      prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
      f1   <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0) 2 * prec * sens / (prec + sens) else NA_real_

      denom_mcc <- sqrt(as.numeric(TP + FP) * as.numeric(TP + FN) * as.numeric(TN + FP) * as.numeric(TN + FN))
      mcc <- if (!is.na(denom_mcc) && denom_mcc > 0) ((TP * TN) - (FP * FN)) / denom_mcc else NA_real_

      diff_vec <- est_vec - truth_vec

      abs_bias <- mean(abs(diff_vec))
      rel_bias <- mean(abs(diff_vec) / (abs(truth_vec) + eps))
      rmse     <- sqrt(mean(diff_vec^2))
      mae      <- mean(abs(diff_vec))
      sd_error <- sqrt(mean((diff_vec - mean(diff_vec))^2))

      if (any(true_nz)) {
        diff_nz       <- diff_vec[true_nz]
        true_nz_vals  <- truth_vec[true_nz]
        abs_bias_nz   <- mean(abs(diff_nz))
        rel_bias_nz   <- mean(abs(diff_nz) / (abs(true_nz_vals) + eps))
        sd_error_nz   <- sqrt(mean((diff_nz - mean(diff_nz))^2))
      } else {
        abs_bias_nz <- NA_real_
        rel_bias_nz <- NA_real_
        sd_error_nz <- NA_real_
      }

      out_list[[length(out_list) + 1]] <- data.frame(
        model        = model_label,
        effect       = "common_tvp",
        subject      = NA_integer_,
        TP           = TP, FP = FP, TN = TN, FN = FN,
        sensitivity  = sens,
        specificity  = spec,
        precision    = prec,
        F1           = f1,
        MCC          = mcc,
        abs_bias     = abs_bias,
        rel_bias     = rel_bias,
        abs_bias_nz  = abs_bias_nz,
        rel_bias_nz  = rel_bias_nz,
        sd_error     = sd_error,
        sd_error_nz  = sd_error_nz,
        RMSE         = rmse,
        MAE          = mae,
        MSFE_min     = msfe_min,
        stringsAsFactors = FALSE
      )
    }
  }

  ## ---- unique_tvp (for TVP models, per subject) ----
  if (is_tvp && length(truth_unique_tvp) > 0 && length(est_unique_tvp) > 0) {
    # Unique TVP: list of list of matrices (subjects x time)
    k_tvp <- min(length(truth_unique_tvp), length(est_unique_tvp))

    for (ii in seq_len(k_tvp)) {
      truth_tvp_list <- truth_unique_tvp[[ii]]
      est_tvp_list   <- est_unique_tvp[[ii]]

      # Only proceed if both are lists
      if (is.list(truth_tvp_list) && is.list(est_tvp_list)) {
        n_t <- min(length(truth_tvp_list), length(est_tvp_list))

        if (n_t > 0) {
          truth_mats <- Filter(is_mat, truth_tvp_list[1:n_t])
          est_mats   <- Filter(is_mat, est_tvp_list[1:n_t])

          if (length(truth_mats) > 0 && length(est_mats) > 0) {
            # Flatten across time
            truth_vec <- as.numeric(do.call(cbind, truth_mats))
            est_vec   <- as.numeric(do.call(cbind, est_mats))

            # Binary edge indicators
            true_nz <- (abs(truth_vec) > 0)
            est_nz  <- (abs(est_vec)  > 0)

            TP <- sum(true_nz & est_nz)
            TN <- sum(!true_nz & !est_nz)
            FP <- sum(!true_nz &  est_nz)
            FN <- sum( true_nz & !est_nz)

            sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
            spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
            prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
            f1   <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0) 2 * prec * sens / (prec + sens) else NA_real_

            denom_mcc <- sqrt(as.numeric(TP + FP) * as.numeric(TP + FN) * as.numeric(TN + FP) * as.numeric(TN + FN))
            mcc <- if (!is.na(denom_mcc) && denom_mcc > 0) ((TP * TN) - (FP * FN)) / denom_mcc else NA_real_

            diff_vec <- est_vec - truth_vec

            abs_bias <- mean(abs(diff_vec))
            rel_bias <- mean(abs(diff_vec) / (abs(truth_vec) + eps))
            rmse     <- sqrt(mean(diff_vec^2))
            mae      <- mean(abs(diff_vec))
            sd_error <- sqrt(mean((diff_vec - mean(diff_vec))^2))

            if (any(true_nz)) {
              diff_nz       <- diff_vec[true_nz]
              true_nz_vals  <- truth_vec[true_nz]
              abs_bias_nz   <- mean(abs(diff_nz))
              rel_bias_nz   <- mean(abs(diff_nz) / (abs(true_nz_vals) + eps))
              sd_error_nz   <- sqrt(mean((diff_nz - mean(diff_nz))^2))
            } else {
              abs_bias_nz <- NA_real_
              rel_bias_nz <- NA_real_
              sd_error_nz <- NA_real_
            }

            out_list[[length(out_list) + 1]] <- data.frame(
              model        = model_label,
              effect       = "tvp",
              subject      = ii,
              TP           = TP, FP = FP, TN = TN, FN = FN,
              sensitivity  = sens,
              specificity  = spec,
              precision    = prec,
              F1           = f1,
              MCC          = mcc,
              abs_bias     = abs_bias,
              rel_bias     = rel_bias,
              abs_bias_nz  = abs_bias_nz,
              rel_bias_nz  = rel_bias_nz,
              sd_error     = sd_error,
              sd_error_nz  = sd_error_nz,
              RMSE         = rmse,
              MAE          = mae,
              MSFE_min     = msfe_min,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  ## ---- intercepts (if requested) ----
  if (intercept) {
    # Common intercept
    if (!is.null(truth_intercepts$common) && !is.null(est_intercepts$common)) {
      out_list[[length(out_list) + 1]] <-
        compute_metrics_for_vector(truth_intercepts$common, est_intercepts$common,
                                   "intercept_common", subject_id = NA_integer_,
                                   eps = eps, msfe_min_global = msfe_min)
    }

    # Unique intercepts (per subject)
    if (length(truth_intercepts$unique) > 0 && length(est_intercepts$unique) > 0) {
      k_ui <- min(length(truth_intercepts$unique), length(est_intercepts$unique))
      for (ii in seq_len(k_ui)) {
        out_list[[length(out_list) + 1]] <-
          compute_metrics_for_vector(truth_intercepts$unique[[ii]], est_intercepts$unique[[ii]],
                                     "intercept_unique", subject_id = ii,
                                     eps = eps, msfe_min_global = msfe_min)
      }
    }

    # Total intercepts (per subject)
    if (length(truth_intercepts$total) > 0 && length(est_intercepts$total) > 0) {
      k_ti <- min(length(truth_intercepts$total), length(est_intercepts$total))
      for (ii in seq_len(k_ti)) {
        out_list[[length(out_list) + 1]] <-
          compute_metrics_for_vector(truth_intercepts$total[[ii]], est_intercepts$total[[ii]],
                                     "intercept_total", subject_id = ii,
                                     eps = eps, msfe_min_global = msfe_min)
      }
    }
  }
  
  ## ---- assemble ----
  if (length(out_list) == 0L) {
    perf_df <- data.frame(
      model = character(), effect = character(), subject = integer(),
      TP = integer(), FP = integer(), TN = integer(), FN = integer(),
      sensitivity = numeric(), specificity = numeric(),
      precision = numeric(), F1 = numeric(), MCC = numeric(),
      abs_bias = numeric(), rel_bias = numeric(),
      abs_bias_nz = numeric(), rel_bias_nz = numeric(),
      sd_error = numeric(), sd_error_nz = numeric(),
      RMSE = numeric(), MAE = numeric(),
      MSFE_min = numeric(),
      stringsAsFactors = FALSE
    )
    return(perf_df)
  }
  
  perf_df <- do.call(rbind, out_list)
  
  ## ---- optional summaries only for present effects ----
  if (any(perf_df$effect == "unique")) {
    perf_df <- rbind(perf_df, add_summary_block(perf_df, "unique", "unique_mean"))
  }
  if (any(perf_df$effect == "subgrp")) {
    perf_df <- rbind(perf_df, add_summary_block(perf_df, "subgrp", "subgrp_mean"))
  }
  if (any(perf_df$effect == "tvp")) {
    perf_df <- rbind(perf_df, add_summary_block(perf_df, "tvp", "tvp_mean"))
  }
  if (any(perf_df$effect == "total")) {
    perf_df <- rbind(perf_df, add_summary_block(perf_df, "total",  "total_mean"))
  }
  if (intercept) {
    if (any(perf_df$effect == "intercept_unique")) {
      perf_df <- rbind(perf_df, add_summary_block(perf_df, "intercept_unique", "intercept_unique_mean"))
    }
    if (any(perf_df$effect == "intercept_total")) {
      perf_df <- rbind(perf_df, add_summary_block(perf_df, "intercept_total", "intercept_total_mean"))
    }
  }
  
  rownames(perf_df) <- NULL
  
  ## ---- keep only averages if requested ----
  if (isTRUE(averages.only)) {
    # Keep summary rows plus 'common', 'common_tvp' and intercept summaries (single aggregates) when present.
    keep_effects <- c("unique_mean", "subgrp_mean", "tvp_mean", "total_mean", "common", "common_tvp",
                      "intercept_common", "intercept_unique_mean", "intercept_total_mean")
    keep_rows <- (perf_df$effect %in% keep_effects) & is.na(perf_df$subject)
    perf_df <- perf_df[keep_rows, , drop = FALSE]
  }
  
  ## ---- reduced output option ----
  if (isTRUE(reduced.output)) {
    # Base columns always present
    keep <- c("model","effect","subject",
              "abs_bias","rel_bias","sd_error","RMSE")

    # Add additional columns if present
    if ("abs_bias_nz" %in% colnames(perf_df)) {
      keep <- c(keep, "abs_bias_nz", "rel_bias_nz", "sd_error_nz")
    }
    if ("sensitivity" %in% colnames(perf_df)) {
      keep <- c(keep, "sensitivity", "specificity", "MCC")
    }

    perf_df <- perf_df[, intersect(keep, colnames(perf_df)), drop = FALSE]
  } else {
    # Ensure 'model' is the first column
    ord <- c("model","effect","subject",
             setdiff(colnames(perf_df), c("model","effect","subject")))
    perf_df <- perf_df[, ord, drop = FALSE]
  }
  
  perf_df
}

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
#' @export
eval_multivar_performance_older <- function(sim_obj,
                                      fit_obj,
                                      eps = 1e-8,
                                      reduced.output = TRUE,
                                      averages.only = TRUE,
                                      label = NULL) {
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
    abs_bias <- mean(abs(diff_vec))
    rel_bias <- mean(abs(diff_vec) / (abs(true_vec) + eps))
    rmse     <- sqrt(mean(diff_vec^2))
    mae      <- mean(abs(diff_vec))
    
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
      RMSE         = rmse,
      MAE          = mae,
      MSFE_min     = msfe_min_global,
      stringsAsFactors = FALSE
    )
  }
  
  add_summary_block <- function(df, which_effect, label) {
    sub <- df[df$effect == which_effect & !is.na(df$subject), , drop = FALSE]
    if (nrow(sub) == 0) return(NULL)
    num_cols <- c("TP","FP","TN","FN",
                  "sensitivity","specificity","precision","F1","MCC",
                  "abs_bias","rel_bias","RMSE","MAE","MSFE_min")
    means <- colMeans(sub[, num_cols, drop = FALSE], na.rm = TRUE)
    data.frame(model = model_label,
               effect = label, subject = NA_integer_, t(means),
               stringsAsFactors = FALSE, row.names = NULL)
  }
  
  ## ---- pull truth (optional) ----
  truth_common  <- if (!is.null(sim_obj$mat_com) && is_mat(sim_obj$mat_com)) sim_obj$mat_com else NULL
  truth_uniqueL <- Filter(is_mat, as_list(sim_obj$mat_ind_unique))
  truth_totalL  <- Filter(is_mat, as_list(sim_obj$mat_ind_final))
  
  ## ---- pull estimates (optional) ----
  est_common  <- if (!is.null(fit_obj$mats$common) && is_mat(fit_obj$mats$common)) fit_obj$mats$common else NULL
  est_uniqueL <- Filter(is_mat, as_list(fit_obj$mats$unique))
  est_totalL  <- Filter(is_mat, as_list(fit_obj$mats$total))
  
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
  
  ## ---- total (match lengths safely) ----
  if (length(truth_totalL) > 0 && length(est_totalL) > 0) {
    k_t <- min(length(truth_totalL), length(est_totalL))
    for (ii in seq_len(k_t)) {
      out_list[[length(out_list) + 1]] <-
        compute_metrics_for_pair(truth_totalL[[ii]], est_totalL[[ii]],
                                 "total", subject_id = ii, eps = eps,
                                 msfe_min_global = msfe_min)
    }
  }
  
  ## ---- assemble ----
  if (length(out_list) == 0L) {
    perf_df <- data.frame(
      model = character(), effect = character(), subject = integer(),
      TP = integer(), FP = integer(), TN = integer(), FN = integer(),
      sensitivity = numeric(), specificity = numeric(),
      precision = numeric(), F1 = numeric(), MCC = numeric(),
      abs_bias = numeric(), rel_bias = numeric(), RMSE = numeric(), MAE = numeric(),
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
  if (any(perf_df$effect == "total")) {
    perf_df <- rbind(perf_df, add_summary_block(perf_df, "total",  "total_mean"))
  }
  
  rownames(perf_df) <- NULL
  
  ## ---- keep only averages if requested ----
  if (isTRUE(averages.only)) {
    # Keep summary rows plus 'common' (single aggregate) when present.
    keep_rows <- (perf_df$effect %in% c("unique_mean","total_mean","common")) & is.na(perf_df$subject)
    perf_df <- perf_df[keep_rows, , drop = FALSE]
  }
  
  ## ---- reduced output option ----
  if (isTRUE(reduced.output)) {
    keep <- c("model","effect","subject","abs_bias","rel_bias","RMSE","sensitivity","specificity","MCC")
    perf_df <- perf_df[, keep, drop = FALSE]
  } else {
    # Ensure 'model' is the first column
    ord <- c("model","effect","subject", setdiff(colnames(perf_df), c("model","effect","subject")))
    perf_df <- perf_df[, ord, drop = FALSE]
  }
  
  perf_df
}

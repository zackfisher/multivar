library("multivar")

set.seed(456)

d <- 4
k <- 1
n <- 240

# Create time-varying dynamics
phi_period1 <- matrix(c(
  0.3,  0.4,  0.0,  0.0,
  0.0,  0.5,  0.3,  0.0,
  0.4,  0.0,  0.2,  0.0,
  0.0,  0.3,  0.0,  0.4
), nrow = d, byrow = TRUE)

phi_period2 <- matrix(c(
  0.4,  0.3,  0.0,  0.2,
  0.0,  0.6,  0.4,  0.0,
  0.3,  0.0,  0.3,  0.0,
  0.0,  0.4,  0.0,  0.5
), nrow = d, byrow = TRUE)

phi_period3 <- matrix(c(
  0.5,  0.2,  0.0,  0.0,
  0.0,  0.4,  0.5,  0.0,
  0.5,  0.0,  0.4,  0.2,
  0.0,  0.5,  0.0,  0.3
), nrow = d, byrow = TRUE)

# Create list of phi matrices for each time point
phi_list <- c(
  rep(list(phi_period1), 80),
  rep(list(phi_period2), 80),
  rep(list(phi_period3), 80)
)

# Generate data
data1 <- multivar:::var_sim_growth(
  n = n,
  phi = phi_list,
  sigma = diag(d),
  intercept = rep(0, d)
)

data <- list(data1)
breaks <- list(c(80, 160))

# Fit the model
object <- multivar::constructModel(
  data,
  intercept = FALSE,
  tvp = TRUE,
  breaks = breaks,
  cv = "blocked",
  nfolds = 5,
  weightest = "lasso",
  lambda_choice = "lambda.min"
)
fit <- multivar::cv.multivar(object)

#-------------------------------------------------------#
context("test06: K=1 TVP runs without error")
#-------------------------------------------------------#

test_that("K=1 TVP model runs successfully", {
  expect_false(is.null(fit))
  expect_false(is.null(fit$mats))
})

test_that("K=1 TVP model returns required components", {
  expect_false(is.null(fit$mats$common))
  expect_false(is.null(fit$mats$unique))
  expect_false(is.null(fit$mats$total))
})

#-------------------------------------------------------#
context("test06: K=1 TVP structure")
#-------------------------------------------------------#

test_that("common is matrix for K=1 TVP (base effects)", {
  expect_true(is.matrix(fit$mats$common))
  expect_equal(dim(fit$mats$common), c(d, d))
})

test_that("unique is list containing base effects matrix", {
  expect_true(is.list(fit$mats$unique))
  expect_true(is.matrix(fit$mats$unique[[1]]))
  expect_equal(dim(fit$mats$unique[[1]]), c(d, d))
})

test_that("tvp is list of lists for TVP effects", {
  expect_true(is.list(fit$mats$tvp))
  expect_true(is.list(fit$mats$tvp[[1]]))
})

test_that("total is list of lists for total effects", {
  expect_true(is.list(fit$mats$total))
  expect_true(is.list(fit$mats$total[[1]]))
})

test_that("number of periods is correct", {
  num_periods <- length(breaks[[1]]) + 1
  expect_equal(length(fit$mats$tvp[[1]]), num_periods)
  expect_equal(length(fit$mats$total[[1]]), num_periods)
})

#-------------------------------------------------------#
context("test06: K=1 common = unique (base effects)")
#-------------------------------------------------------#

test_that("common equals unique for K=1 (base effects)", {
  # For K=1, common and unique[[1]] should be the same (base effects)
  expect_equal(fit$mats$common, fit$mats$unique[[1]])
})

test_that("common is base effects (shared across periods)", {
  # For K=1, common represents the base effects
  # It should be a single d x d matrix
  expect_equal(dim(fit$mats$common), c(d, d))
})

#-------------------------------------------------------#
context("test06: K=1 TVP matrix dimensions")
#-------------------------------------------------------#

test_that("each period matrix has correct dimensions", {
  num_periods <- length(fit$mats$total[[1]])
  for (p in 1:num_periods) {
    expect_equal(dim(fit$mats$total[[1]][[p]]), c(d, d))
  }
})

#-------------------------------------------------------#
context("test06: K=1 TVP no subgroup effects")
#-------------------------------------------------------#

test_that("subgroup effects are NULL for K=1", {
  expect_true(is.null(fit$mats$subgrp))
})

#-------------------------------------------------------#
context("test06: K=1 TVP common effects correct")
#-------------------------------------------------------#

expect_equal_to_reference(
  fit$mats$common, "rds/test06_k1_tvp_common.rds"
)

#-------------------------------------------------------#
context("test06: K=1 TVP unique effects correct")
#-------------------------------------------------------#

expect_equal_to_reference(
  fit$mats$unique, "rds/test06_k1_tvp_unique.rds"
)

#-------------------------------------------------------#
context("test06: K=1 TVP total effects correct")
#-------------------------------------------------------#

expect_equal_to_reference(
  fit$mats$total, "rds/test06_k1_tvp_total.rds"
)

#-------------------------------------------------------#
context("test06: K=1 TVP results are deterministic")
#-------------------------------------------------------#

test_that("results are identical to saved snapshots", {
  expect_identical(fit$mats$common, readRDS("rds/test06_k1_tvp_common.rds"))
  expect_identical(fit$mats$unique, readRDS("rds/test06_k1_tvp_unique.rds"))
  expect_identical(fit$mats$total, readRDS("rds/test06_k1_tvp_total.rds"))
})

#-------------------------------------------------------#
context("test06: K=1 TVP period-specific dynamics")
#-------------------------------------------------------#

test_that("period 1 dynamics are reasonable", {
  period1 <- fit$mats$total[[1]][[1]]
  # Should have some non-zero entries
  expect_gt(sum(abs(period1) > 0), 0)
  # Should be d x d
  expect_equal(dim(period1), c(d, d))
})

test_that("period 2 dynamics are reasonable", {
  period2 <- fit$mats$total[[1]][[2]]
  # Should have some non-zero entries
  expect_gt(sum(abs(period2) > 0), 0)
  # Should be d x d
  expect_equal(dim(period2), c(d, d))
})

test_that("period 3 dynamics are reasonable", {
  period3 <- fit$mats$total[[1]][[3]]
  # Should have some non-zero entries
  expect_gt(sum(abs(period3) > 0), 0)
  # Should be d x d
  expect_equal(dim(period3), c(d, d))
})

test_that("dynamics differ across periods", {
  # Periods should have different dynamics (at least some differences)
  period1 <- fit$mats$total[[1]][[1]]
  period2 <- fit$mats$total[[1]][[2]]
  period3 <- fit$mats$total[[1]][[3]]

  # Not all periods should be identical
  all_same <- identical(period1, period2) && identical(period2, period3)
  expect_false(all_same)
})

#-------------------------------------------------------#
context("test06: K=1 TVP recovery metrics")
#-------------------------------------------------------#

test_that("period 1 has good structure recovery", {
  true_mat <- phi_period1
  est_mat <- fit$mats$total[[1]][[1]]

  true_vec <- as.vector(true_mat)
  est_vec <- as.vector(est_mat)

  true_nz <- abs(true_vec) > 0
  est_nz <- abs(est_vec) > 0

  TP <- sum(true_nz & est_nz)
  TN <- sum(!true_nz & !est_nz)
  FP <- sum(!true_nz & est_nz)
  FN <- sum(true_nz & !est_nz)

  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA

  # Should find most true edges
  expect_gt(sens, 0.8)
  # Should have reasonable specificity
  expect_gt(spec, 0.5)
})

test_that("period 2 has good structure recovery", {
  true_mat <- phi_period2
  est_mat <- fit$mats$total[[1]][[2]]

  true_vec <- as.vector(true_mat)
  est_vec <- as.vector(est_mat)

  true_nz <- abs(true_vec) > 0
  est_nz <- abs(est_vec) > 0

  TP <- sum(true_nz & est_nz)
  TN <- sum(!true_nz & !est_nz)
  FP <- sum(!true_nz & est_nz)
  FN <- sum(true_nz & !est_nz)

  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA

  # Should find most true edges
  expect_gt(sens, 0.8)
  # Should have reasonable specificity
  expect_gt(spec, 0.5)
})

test_that("period 3 has good structure recovery", {
  true_mat <- phi_period3
  est_mat <- fit$mats$total[[1]][[3]]

  true_vec <- as.vector(true_mat)
  est_vec <- as.vector(est_mat)

  true_nz <- abs(true_vec) > 0
  est_nz <- abs(est_vec) > 0

  TP <- sum(true_nz & est_nz)
  TN <- sum(!true_nz & !est_nz)
  FP <- sum(!true_nz & est_nz)
  FN <- sum(true_nz & !est_nz)

  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA

  # Should find most true edges
  expect_gt(sens, 0.8)
  # Should have reasonable specificity
  expect_gt(spec, 0.5)
})

#-------------------------------------------------------#
context("test06: eval_multivar_performance with TVP")
#-------------------------------------------------------#

test_that("eval function works with K=1 TVP", {
  sim_obj <- list(
    k = 1,
    mat_ind_final = list(list(phi_period1, phi_period2, phi_period3))
  )

  result <- multivar::eval_multivar_performance(
    sim_obj = sim_obj,
    fit_obj = fit,
    intercept = FALSE,
    reduced.output = FALSE,
    averages.only = FALSE
  )

  expect_false(is.null(result))
  expect_true(is.data.frame(result))
})

test_that("eval provides period-specific metrics for TVP", {
  sim_obj <- list(
    k = 1,
    mat_ind_final = list(list(phi_period1, phi_period2, phi_period3))
  )

  result <- multivar::eval_multivar_performance(
    sim_obj = sim_obj,
    fit_obj = fit,
    intercept = FALSE,
    reduced.output = FALSE,
    averages.only = FALSE
  )

  # Should have period-specific effects
  expect_true("total_period1" %in% result$effect)
  expect_true("total_period2" %in% result$effect)
  expect_true("total_period3" %in% result$effect)
})

test_that("eval provides overall aggregated metrics for TVP", {
  sim_obj <- list(
    k = 1,
    mat_ind_final = list(list(phi_period1, phi_period2, phi_period3))
  )

  result <- multivar::eval_multivar_performance(
    sim_obj = sim_obj,
    fit_obj = fit,
    intercept = FALSE,
    reduced.output = FALSE,
    averages.only = FALSE
  )

  # Should have overall total effect
  expect_true("total" %in% result$effect)

  # Overall metrics should be aggregated
  overall <- result[result$effect == "total", ]
  expect_equal(nrow(overall), 1)
  expect_false(is.na(overall$MCC))
  expect_false(is.na(overall$sensitivity))
  expect_false(is.na(overall$specificity))
})

test_that("period-specific metrics are reasonable", {
  sim_obj <- list(
    k = 1,
    mat_ind_final = list(list(phi_period1, phi_period2, phi_period3))
  )

  result <- multivar::eval_multivar_performance(
    sim_obj = sim_obj,
    fit_obj = fit,
    intercept = FALSE,
    reduced.output = FALSE,
    averages.only = FALSE
  )

  for (p in 1:3) {
    period_res <- result[result$effect == paste0("total_period", p), ]
    expect_equal(nrow(period_res), 1)

    # Should have good structure recovery
    expect_gt(period_res$MCC, 0.5)
    expect_gt(period_res$sensitivity, 0.8)
    expect_gt(period_res$specificity, 0.5)

    # Should have reasonable error metrics
    expect_lt(period_res$RMSE, 0.2)
  }
})

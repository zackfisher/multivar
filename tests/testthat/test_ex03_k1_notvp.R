library("multivar")

set.seed(123)

d <- 4
k <- 1
n <- 100

# Simulate data for a single subject
sim <- multivar::multivar_sim(
  k = k,
  d = d,
  n = n,
  prop_fill_com = 0.3,
  prop_fill_ind = 0.15,
  unique_overlap = TRUE,
  sigma = diag(d),
  lb = 0.1,
  ub = 0.9
)

# Fit the model
object <- multivar::constructModel(sim$data)
fit <- multivar::cv.multivar(object)

#-------------------------------------------------------#
context("test03: K=1 model runs without error")
#-------------------------------------------------------#

test_that("K=1 model runs successfully", {
  expect_false(is.null(fit))
  expect_false(is.null(fit$mats))
})

test_that("K=1 model returns required components", {
  expect_false(is.null(fit$mats$common))
  expect_false(is.null(fit$mats$unique))
  expect_false(is.null(fit$mats$total))
})

#-------------------------------------------------------#
context("test03: K=1 common = unique = total")
#-------------------------------------------------------#

test_that("common equals total for K=1", {
  expect_equal(fit$mats$common, fit$mats$total[[1]])
})

test_that("unique equals total for K=1", {
  expect_equal(fit$mats$unique[[1]], fit$mats$total[[1]])
})

test_that("common equals unique for K=1", {
  expect_equal(fit$mats$common, fit$mats$unique[[1]])
})

#-------------------------------------------------------#
context("test03: K=1 matrix dimensions")
#-------------------------------------------------------#

test_that("common matrix has correct dimensions", {
  expect_equal(dim(fit$mats$common), c(d, d))
})

test_that("unique matrix has correct dimensions", {
  expect_equal(length(fit$mats$unique), 1)
  expect_equal(dim(fit$mats$unique[[1]]), c(d, d))
})

test_that("total matrix has correct dimensions", {
  expect_equal(length(fit$mats$total), 1)
  expect_equal(dim(fit$mats$total[[1]]), c(d, d))
})

#-------------------------------------------------------#
context("test03: K=1 no subgroup effects")
#-------------------------------------------------------#

test_that("subgroup effects are NULL for K=1", {
  expect_true(is.null(fit$mats$subgrp))
})

#-------------------------------------------------------#
context("test03: K=1 common effects correct")
#-------------------------------------------------------#

expect_equal_to_reference(
  fit$mats$common, "rds/test03_k1_common_effects.rds"
)

#-------------------------------------------------------#
context("test03: K=1 unique effects correct")
#-------------------------------------------------------#

expect_equal_to_reference(
  fit$mats$unique, "rds/test03_k1_unique_effects.rds"
)

#-------------------------------------------------------#
context("test03: K=1 total effects correct")
#-------------------------------------------------------#

expect_equal_to_reference(
  fit$mats$total, "rds/test03_k1_total_effects.rds"
)

#-------------------------------------------------------#
context("test03: K=1 results are deterministic")
#-------------------------------------------------------#

test_that("results are identical to saved snapshots", {
  expect_identical(fit$mats$common, readRDS("rds/test03_k1_common_effects.rds"))
  expect_identical(fit$mats$unique, readRDS("rds/test03_k1_unique_effects.rds"))
  expect_identical(fit$mats$total, readRDS("rds/test03_k1_total_effects.rds"))
})

#-------------------------------------------------------#
context("test03: K=1 evaluation metrics")
#-------------------------------------------------------#

test_that("eval_multivar_performance works for K=1", {
  eval_results <- multivar::eval_multivar_performance(
    sim_obj = sim,
    fit_obj = fit,
    reduced.output = FALSE,
    averages.only = FALSE
  )

  # For K=1, should only have total effects (no common/unique)
  expect_true("total" %in% eval_results$effect)
  expect_false("common" %in% eval_results$effect)
  expect_false("unique" %in% eval_results$effect)
})

test_that("K=1 total effects have good performance", {
  eval_results <- multivar::eval_multivar_performance(
    sim_obj = sim,
    fit_obj = fit,
    reduced.output = FALSE,
    averages.only = FALSE
  )

  total_stats <- eval_results[eval_results$effect == "total", ]

  # Perfect structure recovery expected
  expect_equal(total_stats$MCC[1], 1.0)
  expect_equal(total_stats$sensitivity[1], 1.0)
  expect_equal(total_stats$specificity[1], 1.0)

  # Good parameter recovery expected
  expect_lt(total_stats$RMSE[1], 0.15)
})

library("multivar")

set.seed(789)

d <- 4
k <- 1
n <- 150

# Generate base dynamics
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

# Create non-zero intercepts
set.seed(100)
true_intercept <- runif(d, -0.5, 0.5)

# Add intercepts to the data
data_with_intercept <- sim$data[[1]]
for (i in 1:nrow(data_with_intercept)) {
  data_with_intercept[i, ] <- data_with_intercept[i, ] + true_intercept
}

sim_with_intercept <- sim
sim_with_intercept$data <- list(data_with_intercept)
sim_with_intercept$intercept <- list(true_intercept)

# Fit the model with intercept
object <- multivar::constructModel(
  sim_with_intercept$data,
  intercept = TRUE,
  weightest = "ols",
  nfolds = 10
)
fit <- multivar::cv.multivar(object)

#-------------------------------------------------------#
context("test05: K=1 intercept model runs without error")
#-------------------------------------------------------#

test_that("K=1 intercept model runs successfully", {
  expect_false(is.null(fit))
  expect_false(is.null(fit$mats))
})

test_that("K=1 intercept penalized model returns required components", {
  expect_false(is.null(fit$mats$common))
  expect_false(is.null(fit$mats$unique))
  expect_false(is.null(fit$mats$total))
})

#-------------------------------------------------------#
context("test05: K=1 common = unique = total")
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
context("test05: K=1 matrix dimensions with intercept")
#-------------------------------------------------------#

test_that("common matrix has correct dimensions (d x d, no intercept column)", {
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
context("test05: K=1 intercepts stored separately")
#-------------------------------------------------------#

test_that("intercepts are in separate list", {
  expect_true(!is.null(fit$mats$intercepts))
  expect_true(!is.null(fit$mats$intercepts$intercepts_total))
})

test_that("intercept values are present", {
  intercepts_total <- fit$mats$intercepts$intercepts_total[[1]]
  # Intercepts should be recovered
  expect_equal(length(intercepts_total), d)
})

#-------------------------------------------------------#
context("test05: K=1 no subgroup effects")
#-------------------------------------------------------#

test_that("subgroup effects are NULL for K=1", {
  expect_true(is.null(fit$mats$subgrp))
})

#-------------------------------------------------------#
context("test05: K=1 common effects correct")
#-------------------------------------------------------#

expect_equal_to_reference(
  fit$mats$common, "rds/test05_k1_intercept_penalized_common.rds"
)

#-------------------------------------------------------#
context("test05: K=1 unique effects correct")
#-------------------------------------------------------#

expect_equal_to_reference(
  fit$mats$unique, "rds/test05_k1_intercept_penalized_unique.rds"
)

#-------------------------------------------------------#
context("test05: K=1 total effects correct")
#-------------------------------------------------------#

expect_equal_to_reference(
  fit$mats$total, "rds/test05_k1_intercept_penalized_total.rds"
)

#-------------------------------------------------------#
context("test05: K=1 results are deterministic")
#-------------------------------------------------------#

test_that("results are identical to saved snapshots", {
  expect_identical(fit$mats$common, readRDS("rds/test05_k1_intercept_penalized_common.rds"))
  expect_identical(fit$mats$unique, readRDS("rds/test05_k1_intercept_penalized_unique.rds"))
  expect_identical(fit$mats$total, readRDS("rds/test05_k1_intercept_penalized_total.rds"))
})

#-------------------------------------------------------#
context("test05: K=1 evaluation metrics with intercept")
#-------------------------------------------------------#

test_that("eval_multivar_performance works for K=1 with intercept", {
  eval_results <- multivar::eval_multivar_performance(
    sim_obj = sim_with_intercept,
    fit_obj = fit,
    intercept = TRUE,
    reduced.output = FALSE,
    averages.only = FALSE
  )

  # For K=1, should only have total effects (no common/unique)
  expect_true("total" %in% eval_results$effect)
  expect_false("common" %in% eval_results$effect)
  expect_false("unique" %in% eval_results$effect)

  # Should have intercept evaluation
  expect_true("intercept_total" %in% eval_results$effect)
})

test_that("K=1 total effects have reasonable performance", {
  eval_results <- multivar::eval_multivar_performance(
    sim_obj = sim_with_intercept,
    fit_obj = fit,
    intercept = TRUE,
    reduced.output = FALSE,
    averages.only = FALSE
  )

  total_stats <- eval_results[eval_results$effect == "total", ]

  # Should have good structure recovery
  expect_gt(total_stats$MCC[1], 0.5)
  expect_gt(total_stats$sensitivity[1], 0.8)
  expect_gt(total_stats$specificity[1], 0.5)

  # Should have reasonable parameter recovery
  expect_lt(total_stats$RMSE[1], 0.2)
})

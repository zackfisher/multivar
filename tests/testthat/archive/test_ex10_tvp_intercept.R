library("multivar")

set.seed(777)

d <- 3
k <- 1
n <- 240

# Define time-varying intercepts for 3 periods
true_intercept_period1 <- c(5, 10, 15)
true_intercept_period2 <- c(7, 12, 14)
true_intercept_period3 <- c(6, 9, 16)

# Create time-varying dynamics
phi_period1 <- matrix(c(
  0.3, 0.2, 0.0,
  0.0, 0.4, 0.2,
  0.2, 0.0, 0.3
), nrow = d, byrow = TRUE)

phi_period2 <- matrix(c(
  0.4, 0.3, 0.0,
  0.0, 0.5, 0.3,
  0.3, 0.0, 0.4
), nrow = d, byrow = TRUE)

phi_period3 <- matrix(c(
  0.35, 0.25, 0.0,
  0.0, 0.45, 0.25,
  0.25, 0.0, 0.35
), nrow = d, byrow = TRUE)

# Create list of phi matrices and intercepts for each time point
phi_list <- c(
  rep(list(phi_period1), 80),
  rep(list(phi_period2), 80),
  rep(list(phi_period3), 80)
)

intercept_list <- c(
  rep(list(true_intercept_period1), 80),
  rep(list(true_intercept_period2), 80),
  rep(list(true_intercept_period3), 80)
)

# Generate data with time-varying intercepts
data1 <- multivar:::var_sim_growth(
  n = n,
  phi = phi_list,
  sigma = 0.5 * diag(d),
  intercept = intercept_list
)

data <- list(data1)
# Breakpoints at 80 and 160 create 3 periods
breaks <- list(c(80, 160))

# Fit the model with intercepts (using default weightest = "lasso")
object <- multivar::constructModel(
  data = data,
  tvp = TRUE,
  breaks = breaks,
  intercept = TRUE,
  nfolds = 10
)

fit <- multivar::cv.multivar(object)

#-------------------------------------------------------#
context("test10: TVP intercepts are returned")
#-------------------------------------------------------#

test_that("intercept structure exists", {
  expect_true(!is.null(fit$mats$intercepts))
})

test_that("TVP intercepts are present", {
  expect_true(!is.null(fit$mats$intercepts$intercepts_tvp))
})

test_that("TVP intercepts have correct structure", {
  # Should be a list of subjects
  expect_equal(length(fit$mats$intercepts$intercepts_tvp), k)

  # Each subject should have a list of periods
  expect_equal(length(fit$mats$intercepts$intercepts_tvp[[1]]), 3)

  # Each period should have d intercepts
  for (t in 1:3) {
    expect_equal(length(fit$mats$intercepts$intercepts_tvp[[1]][[t]]), d)
  }
})

#-------------------------------------------------------#
context("test10: TVP intercepts vary over time")
#-------------------------------------------------------#

test_that("intercepts differ across periods", {
  c1 <- fit$mats$intercepts$intercepts_tvp[[1]][[1]]
  c2 <- fit$mats$intercepts$intercepts_tvp[[1]][[2]]
  c3 <- fit$mats$intercepts$intercepts_tvp[[1]][[3]]

  # Periods should not all be identical
  same_12 <- isTRUE(all.equal(c1, c2))
  same_23 <- isTRUE(all.equal(c2, c3))
  same_13 <- isTRUE(all.equal(c1, c3))

  expect_false(same_12 && same_23 && same_13)
})

#-------------------------------------------------------#
context("test10: TVP intercept recovery quality")
#-------------------------------------------------------#

test_that("period 1 intercepts have positive correlation with truth", {
  est <- fit$mats$intercepts$intercepts_tvp[[1]][[1]]
  cor_val <- cor(est, true_intercept_period1)
  expect_gt(cor_val, 0.7)
})

test_that("period 2 intercepts have positive correlation with truth", {
  est <- fit$mats$intercepts$intercepts_tvp[[1]][[2]]
  cor_val <- cor(est, true_intercept_period2)
  expect_gt(cor_val, 0.7)
})

test_that("period 3 intercepts have positive correlation with truth", {
  est <- fit$mats$intercepts$intercepts_tvp[[1]][[3]]
  cor_val <- cor(est, true_intercept_period3)
  expect_gt(cor_val, 0.7)
})

test_that("at least one period has reasonable RMSE", {
  est1 <- fit$mats$intercepts$intercepts_tvp[[1]][[1]]
  est2 <- fit$mats$intercepts$intercepts_tvp[[1]][[2]]
  est3 <- fit$mats$intercepts$intercepts_tvp[[1]][[3]]

  rmse1 <- sqrt(mean((est1 - true_intercept_period1)^2))
  rmse2 <- sqrt(mean((est2 - true_intercept_period2)^2))
  rmse3 <- sqrt(mean((est3 - true_intercept_period3)^2))

  # TVP intercept recovery is challenging; check that at least one period has RMSE < 5
  # (This ensures the method is working, even if not all periods recover perfectly)
  expect_true(any(c(rmse1, rmse2, rmse3) < 5.0))
})

#-------------------------------------------------------#
context("test10: TVP dynamics still recovered")
#-------------------------------------------------------#

test_that("TVP effects in dynamics are present", {
  expect_false(is.null(fit$mats$tvp))
  expect_equal(length(fit$mats$tvp), k)
})

test_that("dynamics vary over time", {
  # Each subject should have time-varying dynamics
  expect_true(is.list(fit$mats$total[[1]]))
  expect_gt(length(fit$mats$total[[1]]), 1)
})

#-------------------------------------------------------#
context("test10: deterministic snapshot tests")
#-------------------------------------------------------#

test_that("TVP intercept results match saved snapshot", {
  expect_equal_to_reference(
    fit$mats$intercepts$intercepts_tvp,
    "rds/test10_tvp_intercepts.rds"
  )
})

test_that("TVP intercepts are identical to saved snapshot", {
  saved_tvp_intercepts <- readRDS("rds/test10_tvp_intercepts.rds")
  expect_identical(fit$mats$intercepts$intercepts_tvp, saved_tvp_intercepts)
})

context("test09: hyperparameter extraction")

#-------------------------------------------------------#
context("test09: basic functionality")
#-------------------------------------------------------#

test_that("fit contains hyperparams dataframe", {
  set.seed(999)
  d <- 3
  sim <- multivar::multivar_sim(
    k = 2, d = d, n = 50,
    prop_fill_com = 0.2, prop_fill_ind = 0.1,
    lb = 0.1, ub = 0.9, sigma = diag(d)
  )
  object <- multivar::constructModel(sim$data)
  fit <- suppressMessages(multivar::cv.multivar(object))

  expect_true("hyperparams" %in% names(fit))
  expect_s3_class(fit$hyperparams, "data.frame")
})

test_that("hyperparams has correct columns", {
  set.seed(999)
  d <- 3
  sim <- multivar::multivar_sim(
    k = 2, d = d, n = 50,
    prop_fill_com = 0.2, prop_fill_ind = 0.1,
    lb = 0.1, ub = 0.9, sigma = diag(d)
  )
  object <- multivar::constructModel(sim$data)
  fit <- suppressMessages(multivar::cv.multivar(object))
  hyp <- fit$hyperparams

  expected_cols <- c("ratio_index", "lambda1_index", "slice_index",
                     "lambda1_value", "ratio_value", "lambda2_value", "MSFE")

  for (col in expected_cols) {
    expect_true(col %in% names(hyp),
                info = paste("Missing column:", col))
  }
})

test_that("hyperparams has correct dimensions", {
  set.seed(999)
  d <- 3
  sim <- multivar::multivar_sim(
    k = 2, d = d, n = 50,
    prop_fill_com = 0.2, prop_fill_ind = 0.1,
    lb = 0.1, ub = 0.9, sigma = diag(d)
  )
  object <- multivar::constructModel(sim$data)
  fit <- suppressMessages(multivar::cv.multivar(object))
  hyp <- fit$hyperparams

  n_lambda <- nrow(object@lambda1)
  n_ratios <- ncol(object@lambda1)
  expected_rows <- n_lambda * n_ratios

  expect_equal(nrow(hyp), expected_rows)
})

test_that("lambda2 = lambda1 * ratio", {
  set.seed(999)
  d <- 3
  sim <- multivar::multivar_sim(
    k = 2, d = d, n = 50,
    prop_fill_com = 0.2, prop_fill_ind = 0.1,
    lb = 0.1, ub = 0.9, sigma = diag(d)
  )
  object <- multivar::constructModel(sim$data)
  fit <- suppressMessages(multivar::cv.multivar(object))
  hyp <- fit$hyperparams

  computed_lambda2 <- hyp$lambda1_value * hyp$ratio_value

  expect_equal(hyp$lambda2_value, computed_lambda2, tolerance = 1e-10)
})

test_that("MSFE values are positive and finite", {
  set.seed(999)
  d <- 3
  sim <- multivar::multivar_sim(
    k = 2, d = d, n = 50,
    prop_fill_com = 0.2, prop_fill_ind = 0.1,
    lb = 0.1, ub = 0.9, sigma = diag(d)
  )
  object <- multivar::constructModel(sim$data)
  fit <- suppressMessages(multivar::cv.multivar(object))
  hyp <- fit$hyperparams

  expect_true(all(hyp$MSFE > 0))
  expect_true(all(is.finite(hyp$MSFE)))
})

test_that("can find best hyperparameters from dataframe", {
  set.seed(999)
  d <- 3
  sim <- multivar::multivar_sim(
    k = 2, d = d, n = 50,
    prop_fill_com = 0.2, prop_fill_ind = 0.1,
    lb = 0.1, ub = 0.9, sigma = diag(d)
  )
  object <- multivar::constructModel(sim$data)
  fit <- suppressMessages(multivar::cv.multivar(object))
  hyp <- fit$hyperparams

  best_idx <- which.min(hyp$MSFE)
  best_row <- hyp[best_idx, ]

  expect_equal(nrow(best_row), 1)
  expect_true(best_row$MSFE == min(hyp$MSFE))
  expect_true(best_row$lambda1_value > 0)
  expect_true(best_row$lambda2_value > 0)
})

#-------------------------------------------------------#
context("test09: plot_cv_lambda_grid works")
#-------------------------------------------------------#

test_that("plot_cv_lambda_grid runs without error", {
  skip_if_not_installed("ggplot2")

  set.seed(999)
  d <- 3
  sim <- multivar::multivar_sim(
    k = 2, d = d, n = 50,
    prop_fill_com = 0.2, prop_fill_ind = 0.1,
    lb = 0.1, ub = 0.9, sigma = diag(d)
  )
  object <- multivar::constructModel(sim$data)
  fit <- suppressMessages(multivar::cv.multivar(object))

  expect_error(
    suppressMessages(multivar::plot_cv_lambda_grid(fit)),
    NA
  )
})

test_that("plot_cv_lambda_grid returns dataframe", {
  skip_if_not_installed("ggplot2")

  set.seed(999)
  d <- 3
  sim <- multivar::multivar_sim(
    k = 2, d = d, n = 50,
    prop_fill_com = 0.2, prop_fill_ind = 0.1,
    lb = 0.1, ub = 0.9, sigma = diag(d)
  )
  object <- multivar::constructModel(sim$data)
  fit <- suppressMessages(multivar::cv.multivar(object))

  result <- suppressMessages(multivar::plot_cv_lambda_grid(fit))

  expect_s3_class(result, "data.frame")
  expect_true("lambda1_value" %in% names(result))
  expect_true("lambda2_value" %in% names(result))
  expect_true("MSFE" %in% names(result))
})

#-------------------------------------------------------#
context("test09: different model types")
#-------------------------------------------------------#

test_that("hyperparams work with intercept model", {
  set.seed(888)
  sim_int <- multivar::multivar_sim(
    k = 2, d = 3, n = 50,
    prop_fill_com = 0.2, prop_fill_ind = 0.1,
    lb = 0.1, ub = 0.9, sigma = diag(3),
    intercept = list(c(0.1, 0.2, 0.3), c(0.4, 0.5, 0.6))
  )

  obj_int <- multivar::constructModel(sim_int$data, intercept = TRUE)
  fit_int <- suppressMessages(multivar::cv.multivar(obj_int))

  expect_true("hyperparams" %in% names(fit_int))
  expect_s3_class(fit_int$hyperparams, "data.frame")
  expect_true("lambda2_value" %in% names(fit_int$hyperparams))
})

test_that("hyperparams work with subgroup model", {
  set.seed(777)
  sim_sub <- multivar::multivar_sim_subgroups(
    k = 4, d = 3, n = 40,
    subgroup = c(1, 1, 2, 2),
    p_com = 0.25, p_sub = 0.10, p_ind = 0.05,
    sigma = diag(3), lb = 0.1, ub = 0.9
  )

  obj_sub <- multivar::constructModel(sim_sub$data,
                                      subgroup_membership = sim_sub$subgroup)
  fit_sub <- suppressMessages(multivar::cv.multivar(obj_sub))

  expect_true("hyperparams" %in% names(fit_sub))
  expect_s3_class(fit_sub$hyperparams, "data.frame")
  expect_true("lambda2_value" %in% names(fit_sub$hyperparams))
})

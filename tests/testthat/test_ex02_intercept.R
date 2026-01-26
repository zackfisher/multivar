library("multivar")

set.seed(123)

d <- 3
k <- 4
n <- 150

# Define true intercepts - different for each individual
true_intercepts <- list(
  c(0.5, -0.3, 0.2),   # Subject 1

c(0.4, -0.2, 0.3),   # Subject 2
  c(0.6, -0.4, 0.1),   # Subject 3
  c(0.3, -0.1, 0.4)    # Subject 4
)

# Simulate data with intercepts
sim <- multivar::multivar_sim(
  k = k,
  d = d,
  n = n,
  prop_fill_com = 0.3,
  prop_fill_ind = 0.1,
  unique_overlap = TRUE,
  sigma = diag(d),
  lb = 0.1,
  ub = 0.5,
  intercept = true_intercepts
)

#-------------------------------------------------------#
context("test02: intercept model runs without error")
#-------------------------------------------------------#

test_that("model with intercept=TRUE runs without error", {
  object <- multivar::constructModel(
    sim$data,
    intercept = TRUE,
    weightest = "ols"
  )
  expect_s4_class(object, "multivar")
  expect_true(object@intercept)
})

test_that("cv.multivar runs with intercept=TRUE", {
  object <- multivar::constructModel(
    sim$data,
    intercept = TRUE,
    weightest = "ols"
  )
  fit <- multivar::cv.multivar(object)
  expect_true(!is.null(fit))
  expect_true(!is.null(fit$mats$common))
  expect_true(!is.null(fit$mats$total))
})

#-------------------------------------------------------#
context("test02: intercept column present in output")
#-------------------------------------------------------#

test_that("common matrix includes intercept column when intercept=TRUE", {
  object <- multivar::constructModel(
    sim$data,
    intercept = TRUE,
    weightest = "ols"
  )
  fit <- multivar::cv.multivar(object)

  # Common matrix should have d+1 columns (intercept + d slopes)
  expect_equal(ncol(fit$mats$common), d + 1)
  expect_equal(nrow(fit$mats$common), d)

  # First column should be named "Intercept"
  expect_equal(colnames(fit$mats$common)[1], "Intercept")
})

test_that("individual total matrices include intercept column", {
  object <- multivar::constructModel(
    sim$data,
    intercept = TRUE,
    weightest = "ols"
  )
  fit <- multivar::cv.multivar(object)

  for (i in seq_along(fit$mats$total)) {
    expect_equal(ncol(fit$mats$total[[i]]), d + 1)
    expect_equal(nrow(fit$mats$total[[i]]), d)
    expect_equal(colnames(fit$mats$total[[i]])[1], "Intercept")
  }
})

#-------------------------------------------------------#
context("test02: intercept=FALSE excludes intercept")
#-------------------------------------------------------#

test_that("common matrix excludes intercept when intercept=FALSE", {
  object <- multivar::constructModel(
    sim$data,
    intercept = FALSE,
    weightest = "ols"
  )
  fit <- multivar::cv.multivar(object)

  # Common matrix should have d columns (no intercept)
  expect_equal(ncol(fit$mats$common), d)
  expect_equal(nrow(fit$mats$common), d)
})

#-------------------------------------------------------#
context("test02: intercept recovery quality")
#-------------------------------------------------------#

test_that("recovered intercepts have strong positive correlation with true values", {
  object <- multivar::constructModel(
    sim$data,
    intercept = TRUE,
    weightest = "ols"
  )
  fit <- multivar::cv.multivar(object)

  # Extract intercept column from each total matrix
  recovered_intercepts <- lapply(fit$mats$total, function(mat) mat[, 1])

  # Compare to true intercepts - should be highly correlated
  true_vec <- unlist(true_intercepts)
  recovered_vec <- unlist(recovered_intercepts)

  cor_val <- cor(true_vec, recovered_vec)

  # Expect strong positive correlation (> 0.7)
  expect_gt(cor_val, 0.7)
})

test_that("mean squared error of recovered intercepts is small", {
  object <- multivar::constructModel(
    sim$data,
    intercept = TRUE,
    weightest = "ols"
  )
  fit <- multivar::cv.multivar(object)

  # Extract intercept column from each total matrix
  recovered_intercepts <- lapply(fit$mats$total, function(mat) mat[, 1])

  # Calculate MSE for each subject
  mse_per_subject <- sapply(seq_along(true_intercepts), function(i) {
    mean((true_intercepts[[i]] - recovered_intercepts[[i]])^2)
  })

  # Overall MSE across all subjects and variables
  overall_mse <- mean(mse_per_subject)

  # MSE should be reasonably small (< 0.1 for intercepts in range [-0.5, 0.5])
  expect_lt(overall_mse, 0.1)
})

test_that("recovered intercepts have correct sign for most entries", {
  object <- multivar::constructModel(
    sim$data,
    intercept = TRUE,
    weightest = "ols"
  )
  fit <- multivar::cv.multivar(object)

  # Extract intercept column from each total matrix
  recovered_intercepts <- lapply(fit$mats$total, function(mat) mat[, 1])

  # Check sign agreement
  sign_matches <- sapply(seq_along(true_intercepts), function(i) {
    true_signs <- sign(true_intercepts[[i]])
    recovered_signs <- sign(recovered_intercepts[[i]])
    mean(true_signs == recovered_signs)
  })

  # At least 70% of signs should match on average
  expect_gt(mean(sign_matches), 0.7)
})

test_that("individual subject intercept recovery is reasonable", {
  object <- multivar::constructModel(
    sim$data,
    intercept = TRUE,
    weightest = "ols"
  )
  fit <- multivar::cv.multivar(object)

  # Extract intercept column from each total matrix
  recovered_intercepts <- lapply(fit$mats$total, function(mat) mat[, 1])

  # Check each subject individually
  for (i in seq_along(true_intercepts)) {
    true_i <- true_intercepts[[i]]
    recovered_i <- recovered_intercepts[[i]]

    # Max absolute error should be reasonable (< 0.5 for intercepts in range [-0.5, 0.5])
    max_abs_error <- max(abs(true_i - recovered_i))
    expect_lt(max_abs_error, 0.5)
  }
})

#-------------------------------------------------------#
context("test02: dynamics recovery with intercept")
#-------------------------------------------------------#

test_that("dynamics are still recovered when intercept=TRUE", {
  object <- multivar::constructModel(
    sim$data,
    intercept = TRUE,
    weightest = "ols"
  )
  fit <- multivar::cv.multivar(object)

  # Extract slopes (columns 2:(d+1)) from total matrices
  recovered_dynamics <- lapply(fit$mats$total, function(mat) mat[, -1])

  # Compare to true dynamics
  for (i in seq_along(recovered_dynamics)) {
    true_mat <- sim$mat_ind_final[[i]]
    recovered_mat <- recovered_dynamics[[i]]

    # Frobenius norm of difference should be reasonable
    frob_diff <- norm(true_mat - recovered_mat, "F")
    frob_true <- norm(true_mat, "F")

    # Relative error should be less than 100% (recovery is meaningful)
    relative_error <- frob_diff / frob_true
    expect_lt(relative_error, 1.0)
  }
})

#-------------------------------------------------------#
context("test02: pen_common_intercept and pen_unique_intercept")
#-------------------------------------------------------#

test_that("pen_common_intercept parameter is accepted", {
  object <- multivar::constructModel(
    sim$data,
    intercept = TRUE,
    pen_common_intercept = FALSE,
    pen_unique_intercept = TRUE,
    weightest = "ols"
  )
  expect_false(object@pen_common_intercept)
  expect_true(object@pen_unique_intercept)
})

test_that("model runs with different intercept penalty settings", {
  # Default: pen_common_intercept=FALSE, pen_unique_intercept=TRUE
  object1 <- multivar::constructModel(
    sim$data,
    intercept = TRUE,
    weightest = "ols"
  )
  fit1 <- multivar::cv.multivar(object1)
  expect_true(!is.null(fit1$mats$common))

  # Both unpenalized
  object2 <- multivar::constructModel(
    sim$data,
    intercept = TRUE,
    pen_common_intercept = FALSE,
    pen_unique_intercept = FALSE,
    weightest = "ols"
  )
  fit2 <- multivar::cv.multivar(object2)
  expect_true(!is.null(fit2$mats$common))
})

#-------------------------------------------------------#
context("test02: intercept results are deterministic")
#-------------------------------------------------------#

test_that("results are identical to saved snapshots", {
  object <- multivar::constructModel(
    sim$data,
    intercept = TRUE,
    weightest = "ols"
  )
  fit <- multivar::cv.multivar(object)

  expect_identical(fit$mats$common, readRDS("rds/test02_intercept_common.rds"))
  expect_identical(fit$mats$unique, readRDS("rds/test02_intercept_unique.rds"))
  expect_identical(fit$mats$total, readRDS("rds/test02_intercept_total.rds"))
})

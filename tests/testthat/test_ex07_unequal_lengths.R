library("multivar")

#-------------------------------------------------------#
# Test: Unequal Time Series Lengths
#-------------------------------------------------------#
# This test verifies that multivar can handle subjects
# with different time series lengths (varying n)
#-------------------------------------------------------#

context("test07: unequal time series lengths")

set.seed(456)

# Create 3 subjects with different time series lengths
# Subject 1: n=80, Subject 2: n=100, Subject 3: n=120
sim1 <- multivar_sim(
  k = 1, d = 4, n = 80,
  prop_fill_com = 0.3,
  prop_fill_ind = 0.15,
  sigma = diag(4),
  lb = 0.1, ub = 0.9
)

sim2 <- multivar_sim(
  k = 1, d = 4, n = 100,
  prop_fill_com = 0.3,
  prop_fill_ind = 0.15,
  sigma = diag(4),
  lb = 0.1, ub = 0.9
)

sim3 <- multivar_sim(
  k = 1, d = 4, n = 120,
  prop_fill_com = 0.3,
  prop_fill_ind = 0.15,
  sigma = diag(4),
  lb = 0.1, ub = 0.9
)

# Combine into a single list (K=3 subjects with varying n)
data_list <- list(
  sim1$data[[1]],
  sim2$data[[1]],
  sim3$data[[1]]
)

#-------------------------------------------------------#
# Test 1: Data dimensions are correct
#-------------------------------------------------------#

test_that("data has correct dimensions for each subject", {
  expect_equal(nrow(data_list[[1]]), 80)
  expect_equal(nrow(data_list[[2]]), 100)
  expect_equal(nrow(data_list[[3]]), 120)
  expect_equal(ncol(data_list[[1]]), 4)
  expect_equal(ncol(data_list[[2]]), 4)
  expect_equal(ncol(data_list[[3]]), 4)
})

#-------------------------------------------------------#
# Test 2: Model construction works with unequal lengths
#-------------------------------------------------------#

test_that("constructModel works with unequal time series lengths", {
  # Should construct without warnings now (common_tvp_effects defaults intelligently)
  expect_silent({
    object <- constructModel(data_list, intercept = FALSE)
  })

  # Verify the object has correct k, d, and n
  object <- constructModel(data_list, intercept = FALSE)
  expect_equal(object@k, 3)

  # d can be a vector of length k (one per subject)
  expect_true(all(object@d == 4))

  expect_equal(length(object@n), 3)
  expect_equal(as.numeric(object@n), c(79, 99, 119))  # n-1 for VAR(1)
})

#-------------------------------------------------------#
# Test 3: Model fitting works
#-------------------------------------------------------#

test_that("cv.multivar works with unequal time series lengths", {
  object <- constructModel(data_list, intercept = FALSE)

  # cv.multivar produces output (progress bar), so don't expect silence
  fit <- suppressMessages(cv.multivar(object))

  # Check that fit object has expected structure
  expect_true(!is.null(fit$mats))
  expect_true(!is.null(fit$mats$common))
  expect_true(!is.null(fit$mats$unique))
  expect_true(!is.null(fit$mats$total))

  # Check that we have results for all 3 subjects
  expect_equal(length(fit$mats$total), 3)

  # Check that all matrices have correct dimensions (4x4)
  expect_equal(dim(fit$mats$common), c(4, 4))
  expect_equal(dim(fit$mats$total[[1]]), c(4, 4))
  expect_equal(dim(fit$mats$total[[2]]), c(4, 4))
  expect_equal(dim(fit$mats$total[[3]]), c(4, 4))
})

#-------------------------------------------------------#
# Test 4: Summary function handles unequal lengths
#-------------------------------------------------------#

test_that("summary_multivar correctly reports unequal lengths", {
  object <- constructModel(data_list, intercept = FALSE)
  fit <- suppressMessages(cv.multivar(object))

  # Capture summary output
  summary_output <- capture.output(summary_multivar(fit))
  summary_text <- paste(summary_output, collapse = "\n")

  # Check that summary mentions varying lengths
  expect_true(grepl("79, 99, 119", summary_text))
  expect_true(grepl("varies by subject", summary_text, ignore.case = TRUE))
  expect_true(grepl("Range:", summary_text))
  expect_true(grepl("Mean:", summary_text))
})

#-------------------------------------------------------#
# Test 5: Print dynamics works
#-------------------------------------------------------#

test_that("print_dynamics works with unequal lengths", {
  object <- constructModel(data_list, intercept = FALSE)
  fit <- suppressMessages(cv.multivar(object))

  # Should not throw an error
  expect_silent({
    capture.output(print_dynamics(fit))
  })
})

#-------------------------------------------------------#
# Test 6: Edge prevalence works
#-------------------------------------------------------#

test_that("print_edge_prevalence works with unequal lengths", {
  object <- constructModel(data_list, intercept = FALSE)
  fit <- suppressMessages(cv.multivar(object))

  # Should not throw an error
  expect_silent({
    capture.output(print_edge_prevalence(fit, type = "proportion"))
  })

  expect_silent({
    capture.output(print_edge_prevalence(fit, type = "count"))
  })
})

#-------------------------------------------------------#
# Test 7: Verify results are stable across runs
#-------------------------------------------------------#

context("test07: unequal time series lengths results are deterministic")

test_that("results are identical to saved snapshots", {
  object <- constructModel(data_list, intercept = FALSE)
  fit <- suppressMessages(cv.multivar(object))

  expect_identical(fit$mats$common, readRDS("rds/test07_unequal_common.rds"))
  expect_identical(fit$mats$total, readRDS("rds/test07_unequal_total.rds"))
})

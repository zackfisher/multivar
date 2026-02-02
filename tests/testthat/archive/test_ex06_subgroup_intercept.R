library("multivar")

set.seed(999)

d <- 3
k_per_group <- 2
n <- 120

# Define true intercepts with subgroup structure
# Common intercept (shared by all)
true_intercept_common <- c(5, 10, 15)

# Subgroup-specific intercepts (shared within subgroup)
true_intercept_subgrp1 <- c(2, -1, 3)
true_intercept_subgrp2 <- c(-2, 1, -3)

# Individual-specific intercepts (unique to each subject)
true_intercepts_unique <- list(
  c(0.5, -0.2, 0.3),   # Subject 1 (subgroup 1)
  c(-0.5, 0.2, -0.3),  # Subject 2 (subgroup 1)
  c(0.3, 0.1, -0.1),   # Subject 3 (subgroup 2)
  c(-0.3, -0.1, 0.1)   # Subject 4 (subgroup 2)
)

# Total intercepts for each subject
true_intercepts_total <- list(
  true_intercept_common + true_intercept_subgrp1 + true_intercepts_unique[[1]],
  true_intercept_common + true_intercept_subgrp1 + true_intercepts_unique[[2]],
  true_intercept_common + true_intercept_subgrp2 + true_intercepts_unique[[3]],
  true_intercept_common + true_intercept_subgrp2 + true_intercepts_unique[[4]]
)

# Simulate data for subgroup 1 (subjects 1-2)
sim1 <- multivar::multivar_sim(
  k = k_per_group,
  d = d,
  n = n,
  prop_fill_com = .25,
  prop_fill_ind = .15,
  unique_overlap = TRUE,
  sigma = 0.5 * diag(d),
  lb = .1,
  ub = .4,
  intercept = list(true_intercepts_total[[1]], true_intercepts_total[[2]])
)

# Simulate data for subgroup 2 (subjects 3-4) with different dynamics
set.seed(1000)
sim2 <- multivar::multivar_sim(
  k = k_per_group,
  d = d,
  n = n,
  prop_fill_com = .25,
  prop_fill_ind = .15,
  unique_overlap = TRUE,
  sigma = 0.5 * diag(d),
  lb = .1,
  ub = .4,
  intercept = list(true_intercepts_total[[3]], true_intercepts_total[[4]])
)

# Combine the two simulations
combined_data <- c(sim1$data, sim2$data)
k <- length(combined_data)

# Subgroup membership: subjects 1-2 in group 1, subjects 3-4 in group 2
subgroup_membership <- c(rep(1, k_per_group), rep(2, k_per_group))

# Fit the model with intercepts
object <- multivar::constructModel(
  combined_data,
  subgroup_membership = subgroup_membership,
  intercept = TRUE,
  weightest = "ols",
  nfolds = 10
)
fit <- multivar::cv.multivar(object)

#-------------------------------------------------------#
context("test06: subgroup intercepts are returned")
#-------------------------------------------------------#

test_that("intercept structure includes subgroup intercepts", {
  expect_true(!is.null(fit$mats$intercepts))
  expect_true(!is.null(fit$mats$intercepts$intercepts_subgrp))
})

test_that("subgroup intercepts have correct length", {
  expect_equal(length(fit$mats$intercepts$intercepts_subgrp), k)
})

test_that("each subgroup intercept is a vector of length d", {
  for (i in seq_len(k)) {
    expect_equal(length(fit$mats$intercepts$intercepts_subgrp[[i]]), d)
  }
})

#-------------------------------------------------------#
context("test06: subjects in same subgroup share intercepts")
#-------------------------------------------------------#

test_that("subjects 1-2 (subgroup 1) have identical subgroup intercepts", {
  expect_equal(
    fit$mats$intercepts$intercepts_subgrp[[1]],
    fit$mats$intercepts$intercepts_subgrp[[2]]
  )
})

test_that("subjects 3-4 (subgroup 2) have identical subgroup intercepts", {
  expect_equal(
    fit$mats$intercepts$intercepts_subgrp[[3]],
    fit$mats$intercepts$intercepts_subgrp[[4]]
  )
})

test_that("subgroup 1 and subgroup 2 have different intercepts", {
  expect_false(isTRUE(all.equal(
    fit$mats$intercepts$intercepts_subgrp[[1]],
    fit$mats$intercepts$intercepts_subgrp[[3]]
  )))
})

#-------------------------------------------------------#
context("test06: intercept decomposition is correct")
#-------------------------------------------------------#

test_that("total = common + subgroup + unique for each subject", {
  for (i in seq_len(k)) {
    reconstructed <- fit$mats$intercepts$intercept_common +
                     fit$mats$intercepts$intercepts_subgrp[[i]] +
                     fit$mats$intercepts$intercepts_unique[[i]]

    expect_equal(
      reconstructed,
      fit$mats$intercepts$intercepts_total[[i]],
      tolerance = 1e-10
    )
  }
})

test_that("sum of unique intercepts is approximately zero", {
  sum_unique <- Reduce("+", fit$mats$intercepts$intercepts_unique)
  expect_true(all(abs(sum_unique) < 1e-10))
})

test_that("common intercept is mean of total intercepts", {
  mean_total <- Reduce("+", fit$mats$intercepts$intercepts_total) / k
  expect_equal(
    fit$mats$intercepts$intercept_common,
    mean_total,
    tolerance = 1e-10
  )
})

#-------------------------------------------------------#
context("test06: intercept recovery quality")
#-------------------------------------------------------#

test_that("common intercept is reasonably close to true common", {
  error <- abs(fit$mats$intercepts$intercept_common - true_intercept_common)
  # With n=120, expect errors < 2
  expect_true(all(error < 2))
})

test_that("subgroup 1 intercepts have positive correlation with truth", {
  est_subgrp1 <- fit$mats$intercepts$intercepts_subgrp[[1]]
  cor_val <- cor(est_subgrp1, true_intercept_subgrp1)
  expect_gt(cor_val, 0.5)
})

test_that("subgroup 2 intercepts have positive correlation with truth", {
  est_subgrp2 <- fit$mats$intercepts$intercepts_subgrp[[3]]
  cor_val <- cor(est_subgrp2, true_intercept_subgrp2)
  expect_gt(cor_val, 0.5)
})

test_that("total intercepts have high correlation with true values", {
  # Flatten to vectors for correlation
  true_vec <- unlist(true_intercepts_total)
  est_vec <- unlist(fit$mats$intercepts$intercepts_total)

  cor_val <- cor(true_vec, est_vec)
  expect_gt(cor_val, 0.8)
})

test_that("RMSE of total intercepts is reasonable", {
  # Compute RMSE for each subject
  rmse_per_subject <- sapply(seq_len(k), function(i) {
    sqrt(mean((true_intercepts_total[[i]] - fit$mats$intercepts$intercepts_total[[i]])^2))
  })

  # With n=120, sigma=0.5, expect RMSE < 2.0
  expect_true(all(rmse_per_subject < 2.0))
})

#-------------------------------------------------------#
context("test06: subgroup dynamics still recovered")
#-------------------------------------------------------#

test_that("subgroup effects in dynamics are present", {
  expect_false(is.null(fit$mats$subgrp))
  expect_equal(length(fit$mats$subgrp), k)
})

test_that("subjects in same subgroup have identical dynamic subgroup effects", {
  # Subgroup 1
  expect_equal(fit$mats$subgrp[[1]], fit$mats$subgrp[[2]])
  # Subgroup 2
  expect_equal(fit$mats$subgrp[[3]], fit$mats$subgrp[[4]])
})

#-------------------------------------------------------#
context("test06: deterministic snapshot tests")
#-------------------------------------------------------#

# Generate snapshots with: saveRDS(fit$mats$intercepts, "rds/test06_subgroup_intercepts.rds")

test_that("intercept results match saved snapshot", {
  expect_equal_to_reference(
    fit$mats$intercepts,
    "rds/test06_subgroup_intercepts.rds"
  )
})

test_that("intercepts are identical to saved snapshot", {
  saved_intercepts <- readRDS("rds/test06_subgroup_intercepts.rds")
  expect_identical(fit$mats$intercepts, saved_intercepts)
})

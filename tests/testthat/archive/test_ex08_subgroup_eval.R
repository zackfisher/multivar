context("test08: subgroup evaluation")

#-------------------------------------------------------#
# Setup: Simulate subgroup data with known ground truth
#-------------------------------------------------------#

set.seed(8080)

# Simulate 6 subjects in 2 subgroups with known structure
sim <- multivar::multivar_sim_subgroups(
  k = 6,
  d = 4,
  n = 80,
  subgroup = c(1, 1, 1, 2, 2, 2),
  p_com = 0.25,
  p_sub = 0.10,
  p_ind = 0.05,
  sigma = diag(4),
  lb = 0.1,
  ub = 0.9,
  intercept = NULL
)

# Fit model with subgroup membership
object <- multivar::constructModel(
  sim$data,
  subgroup_membership = sim$subgroup
)

fit <- suppressMessages(multivar::cv.multivar(object))

#-------------------------------------------------------#
context("test08: subgroup effects are evaluated")
#-------------------------------------------------------#

test_that("eval_multivar_performance works with subgroup data", {
  # Should run without error
  perf <- multivar::eval_multivar_performance(sim, fit)

  expect_s3_class(perf, "data.frame")
  expect_true(nrow(perf) > 0)
})

test_that("subgroup effects appear in performance output", {
  perf <- multivar::eval_multivar_performance(sim, fit, averages.only = FALSE)

  # Check that subgrp effect is present
  expect_true("subgrp" %in% perf$effect)

  # Check that subgrp_mean summary is present
  expect_true("subgrp_mean" %in% perf$effect)
})

test_that("subgrp_mean has valid metrics", {
  perf <- multivar::eval_multivar_performance(sim, fit)

  subgrp_row <- perf[perf$effect == "subgrp_mean", ]

  # Should have exactly one row
  expect_equal(nrow(subgrp_row), 1)

  # Should have NA subject (it's a summary)
  expect_true(is.na(subgrp_row$subject))

  # Should have valid metrics (not all NA)
  expect_false(all(is.na(subgrp_row$sensitivity)))
  expect_false(all(is.na(subgrp_row$specificity)))
  expect_false(all(is.na(subgrp_row$MCC)))
  expect_false(all(is.na(subgrp_row$abs_bias)))
})

test_that("subgroup effects are evaluated per subject", {
  perf <- multivar::eval_multivar_performance(sim, fit, averages.only = FALSE)

  subgrp_rows <- perf[perf$effect == "subgrp" & !is.na(perf$subject), ]

  # Should have one row per subject (6 subjects)
  expect_equal(nrow(subgrp_rows), 6)

  # Should have subjects 1-6
  expect_equal(sort(subgrp_rows$subject), 1:6)
})

test_that("averages.only includes subgrp_mean", {
  perf <- multivar::eval_multivar_performance(sim, fit, averages.only = TRUE)

  # Should still include subgrp_mean
  expect_true("subgrp_mean" %in% perf$effect)

  # Should NOT include individual subgrp rows
  subgrp_individual <- perf[perf$effect == "subgrp" & !is.na(perf$subject), ]
  expect_equal(nrow(subgrp_individual), 0)
})

test_that("all expected effects are present", {
  perf <- multivar::eval_multivar_performance(sim, fit, averages.only = FALSE)

  # Should have: common, unique (6x), subgrp (6x), total (6x),
  # plus summaries: unique_mean, subgrp_mean, total_mean
  expected_effects <- c("common", "unique", "subgrp", "total",
                        "unique_mean", "subgrp_mean", "total_mean")

  for (eff in expected_effects) {
    expect_true(eff %in% perf$effect,
                info = sprintf("Effect '%s' should be present", eff))
  }
})

#-------------------------------------------------------#
context("test08: subgroup structure is correctly simulated")
#-------------------------------------------------------#

test_that("simulation returns all required components", {
  expect_true(!is.null(sim$data))
  expect_true(!is.null(sim$subgroup))
  expect_true(!is.null(sim$mat_com))
  expect_true(!is.null(sim$mat_sub_unique))
  expect_true(!is.null(sim$mat_ind_unique))
  expect_true(!is.null(sim$mat_ind_final))
})

test_that("subgroup membership is correct", {
  expect_equal(length(sim$subgroup), 6)
  expect_equal(sim$subgroup, c(1, 1, 1, 2, 2, 2))
})

test_that("subjects in same subgroup have same subgroup effect", {
  # Subjects 1-3 should have identical subgroup matrices
  expect_equal(sim$mat_sub_unique[[1]], sim$mat_sub_unique[[2]])
  expect_equal(sim$mat_sub_unique[[2]], sim$mat_sub_unique[[3]])

  # Subjects 4-6 should have identical subgroup matrices
  expect_equal(sim$mat_sub_unique[[4]], sim$mat_sub_unique[[5]])
  expect_equal(sim$mat_sub_unique[[5]], sim$mat_sub_unique[[6]])

  # But group 1 and group 2 should differ
  expect_false(identical(sim$mat_sub_unique[[1]], sim$mat_sub_unique[[4]]))
})

test_that("decomposition holds: total = common + subgrp + unique", {
  for (i in 1:6) {
    reconstructed <- sim$mat_com + sim$mat_sub_unique[[i]] + sim$mat_ind_unique[[i]]
    expect_equal(reconstructed, sim$mat_ind_final[[i]], tolerance = 1e-10)
  }
})

test_that("non-overlapping positions: common, subgrp, unique are disjoint", {
  # Check that nonzero positions don't overlap
  for (i in 1:6) {
    com_nz <- which(sim$mat_com != 0)
    sub_nz <- which(sim$mat_sub_unique[[i]] != 0)
    ind_nz <- which(sim$mat_ind_unique[[i]] != 0)

    # Common and subgroup should not overlap
    expect_equal(length(intersect(com_nz, sub_nz)), 0)

    # Common and unique should not overlap
    expect_equal(length(intersect(com_nz, ind_nz)), 0)

    # Subgroup and unique should not overlap
    expect_equal(length(intersect(sub_nz, ind_nz)), 0)
  }
})

#-------------------------------------------------------#
context("test08: results are deterministic")
#-------------------------------------------------------#

test_that("results are identical to saved snapshots", {
  object <- multivar::constructModel(sim$data, subgroup_membership = sim$subgroup, nfolds = 10)
  fit <- suppressMessages(multivar::cv.multivar(object))

  expect_identical(fit$mats$common, readRDS("rds/test08_subgroup_eval_common.rds"))
  expect_identical(fit$mats$subgrp, readRDS("rds/test08_subgroup_eval_subgrp.rds"))
  expect_identical(fit$mats$unique, readRDS("rds/test08_subgroup_eval_unique.rds"))
  expect_identical(fit$mats$total, readRDS("rds/test08_subgroup_eval_total.rds"))
})

library("multivar")

set.seed(456)

d <- 4
k_per_group <- 3
n <- 80

# Simulate data for subgroup 1 (subjects 1-3)
sim1 <- multivar::multivar_sim(
  k = k_per_group,
  d = d,
  n = n,
  prop_fill_com = .2,
  prop_fill_ind = .1,
  unique_overlap = TRUE,
  sigma = diag(d),
  lb = .1,
  ub = .9
)

# Simulate data for subgroup 2 (subjects 4-6) with different seed for different structure
set.seed(789)
sim2 <- multivar::multivar_sim(
  k = k_per_group,
  d = d,
  n = n,
  prop_fill_com = .2,
  prop_fill_ind = .1,
  unique_overlap = TRUE,
  sigma = diag(d),
  lb = .1,
  ub = .9
)

# Combine the two simulations - each subgroup has its own shared structure
combined_data <- c(sim1$data, sim2$data)
k <- length(combined_data)

# Subgroup membership: subjects 1-3 in group 1, subjects 4-6 in group 2
subgroup_membership <- c(rep(1, k_per_group), rep(2, k_per_group))

object <- multivar::constructModel(combined_data, subgroup_membership = subgroup_membership, nfolds = 10)
fit <- multivar::cv.multivar(object)

#-------------------------------------------------------#
context("test01: subgrouping returns subgroup effects")
#-------------------------------------------------------#

test_that("subgroup effects are returned", {
  expect_false(is.null(fit$mats$subgrp))
  expect_equal(length(fit$mats$subgrp), k)
})

test_that("subgroup effects have correct dimensions", {
  for (i in seq_along(fit$mats$subgrp)) {
    expect_equal(dim(fit$mats$subgrp[[i]]), c(d, d))
  }
})

test_that("subjects in same subgroup have same subgroup effects", {
  # Subjects 1-3 should have identical subgroup effect matrices
  expect_equal(fit$mats$subgrp[[1]], fit$mats$subgrp[[2]])
  expect_equal(fit$mats$subgrp[[2]], fit$mats$subgrp[[3]])

  # Subjects 4-6 should have identical subgroup effect matrices
  expect_equal(fit$mats$subgrp[[4]], fit$mats$subgrp[[5]])
  expect_equal(fit$mats$subgrp[[5]], fit$mats$subgrp[[6]])
})

#-------------------------------------------------------#
context("test01: subgrouping common effects correct")
#-------------------------------------------------------#

expect_equal_to_reference(
  fit$mats$common, "rds/test01_common_effects.rds"
)

#-------------------------------------------------------#
context("test01: subgrouping unique effects correct")
#-------------------------------------------------------#

expect_equal_to_reference(
  fit$mats$unique, "rds/test01_unique_effects.rds"
)

#-------------------------------------------------------#
context("test01: subgrouping subgroup effects correct")
#-------------------------------------------------------#

expect_equal_to_reference(
  fit$mats$subgrp, "rds/test01_subgrp_effects.rds"
)

#-------------------------------------------------------#
context("test01: subgrouping total effects correct")
#-------------------------------------------------------#

expect_equal_to_reference(
  fit$mats$total, "rds/test01_total_effects.rds"
)

#-------------------------------------------------------#
context("test01: subgrouping results are deterministic")
#-------------------------------------------------------#

test_that("results are identical to saved snapshots", {
  expect_identical(fit$mats$common, readRDS("rds/test01_common_effects.rds"))
  expect_identical(fit$mats$unique, readRDS("rds/test01_unique_effects.rds"))
  expect_identical(fit$mats$subgrp, readRDS("rds/test01_subgrp_effects.rds"))
  expect_identical(fit$mats$total, readRDS("rds/test01_total_effects.rds"))
})

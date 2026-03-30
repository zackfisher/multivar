library("multivar")

set.seed(123)

d <- 5
k <- 3
n <- 100

sim <- multivar_sim(
  k = k,
  d = d,
  n = n,
  prop_fill_com = .2,
  prop_fill_ind = .1,
  unique_overlap = TRUE,
  sigma = diag(d),
  lb = .1,
  ub = .9
)

#-------------------------------------------------------#
# Step 1: Fit model with CV to find optimal hyperparameters
#-------------------------------------------------------#

object_cv <- constructModel(
  sim$data,
  lambda_choice = "lambda.1se",
  nfolds = 10,
  depth = 1000
)
fit_cv <- cv.multivar(object_cv)

#-------------------------------------------------------#
# Step 2: Extract optimal lambda1 and ratio from fit results
#-------------------------------------------------------#

# Use hyperparams from fit result (has correct lambda1 values)
hyp <- fit_cv$hyperparams
best_row <- hyp[which.min(hyp$MSFE), ]

optimal_lambda1 <- best_row$lambda1_value
optimal_ratio <- best_row$ratio_value

cat("Optimal lambda1:", optimal_lambda1, "\n")
cat("Optimal ratio:", optimal_ratio, "\n")

#-------------------------------------------------------#
# Step 3: Refit model with prespecified lambda1 and ratio
#-------------------------------------------------------#

object_final <- constructModel(
  sim$data,
  lambda1 = optimal_lambda1,
  ratios_unique = optimal_ratio,
  lambda_choice = "lambda.1se",  # Must match original
  nfolds = 10
)
fit_final <- cv.multivar(object_final)

# Extract the final coefficient matrices
fit <- fit_final

#-------------------------------------------------------#
context("test00: common effects correct")
#-------------------------------------------------------#

test_that("common effects match reference", {
  expect_equal_to_reference(
    fit$mats$common, "rds/test00_common_effects.rds"
  )
})

#-------------------------------------------------------#
context("test00: unique effects correct")
#-------------------------------------------------------#

test_that("unique effects match reference", {
  expect_equal_to_reference(
    fit$mats$unique, "rds/test00_unique_effects.rds"
  )
})

#-------------------------------------------------------#
context("test00: total effects correct")
#-------------------------------------------------------#

test_that("total effects match reference", {
  expect_equal_to_reference(
    fit$mats$total, "rds/test00_total_effects.rds"
  )
})

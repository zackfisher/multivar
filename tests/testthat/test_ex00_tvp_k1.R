library("multivar")

d <- 5
T_per_period <- c(200, 200)

phi1 <- matrix(0, d, d)
diag(phi1) <- 0.5
phi1[1, 2] <- 0.3

phi2 <- matrix(0, d, d)
diag(phi2) <- 0.5
phi2[2, 3] <- 0.3

set.seed(123)
sim <- multivar:::var_sim_tvp(
  n = sum(T_per_period),
  sigma = diag(d) * 0.5,
  phi_list = list(phi1, phi2),
  T_per_period = T_per_period,
  type = "piecewise"
)

#-------------------------------------------------------#
# K=1 TVP with common_effects=TRUE, lambda.min
#-------------------------------------------------------#

obj_true_min <- multivar::constructModel(
  sim$data,
  tvp = TRUE,
  breaks = list(sim$breaks),
  lambda_choice = "lambda.min",
  weightest = "lasso",
  common_effects = TRUE
)
fit_true_min <- multivar::cv.multivar(obj_true_min)

context("test00_tvp_k1: common_effects=TRUE, lambda.min - common")
expect_equal_to_reference(
  fit_true_min$mats$common, "rds/test00_tvp_k1_common_TRUE_min_common.rds"
)

context("test00_tvp_k1: common_effects=TRUE, lambda.min - tvp")
expect_equal_to_reference(
  fit_true_min$mats$tvp, "rds/test00_tvp_k1_common_TRUE_min_tvp.rds"
)

#-------------------------------------------------------#
# K=1 TVP with common_effects=TRUE, lambda.1se
#-------------------------------------------------------#

obj_true_1se <- multivar::constructModel(
  sim$data,
  tvp = TRUE,
  breaks = list(sim$breaks),
  lambda_choice = "lambda.1se",
  weightest = "lasso",
  common_effects = TRUE
)
fit_true_1se <- multivar::cv.multivar(obj_true_1se)

context("test00_tvp_k1: common_effects=TRUE, lambda.1se - common")
expect_equal_to_reference(
  fit_true_1se$mats$common, "rds/test00_tvp_k1_common_TRUE_1se_common.rds"
)

context("test00_tvp_k1: common_effects=TRUE, lambda.1se - tvp")
expect_equal_to_reference(
  fit_true_1se$mats$tvp, "rds/test00_tvp_k1_common_TRUE_1se_tvp.rds"
)

#-------------------------------------------------------#
# K=1 TVP with common_effects=FALSE, lambda.min
#-------------------------------------------------------#

obj_false_min <- multivar::constructModel(
  sim$data,
  tvp = TRUE,
  breaks = list(sim$breaks),
  lambda_choice = "lambda.min",
  weightest = "lasso",
  common_effects = FALSE
)
fit_false_min <- multivar::cv.multivar(obj_false_min)

context("test00_tvp_k1: common_effects=FALSE, lambda.min - tvp")
expect_equal_to_reference(
  fit_false_min$mats$tvp, "rds/test00_tvp_k1_common_FALSE_min_tvp.rds"
)

#-------------------------------------------------------#
# K=1 TVP with common_effects=FALSE, lambda.1se
#-------------------------------------------------------#

obj_false_1se <- multivar::constructModel(
  sim$data,
  tvp = TRUE,
  breaks = list(sim$breaks),
  lambda_choice = "lambda.1se",
  weightest = "lasso",
  common_effects = FALSE
)
fit_false_1se <- multivar::cv.multivar(obj_false_1se)

context("test00_tvp_k1: common_effects=FALSE, lambda.1se - tvp")
expect_equal_to_reference(
  fit_false_1se$mats$tvp, "rds/test00_tvp_k1_common_FALSE_1se_tvp.rds"
)

# Test: K=1 TVP with common_effects=FALSE should give nearly identical
# transition matrices as fitting separate period-specific models
# when using the same lambda value.
#
# Note: Small numerical differences are expected due to:
# - Different step sizes in FISTA (eigenvalue of different-sized matrices)
# - Sparse vs dense matrix operations
# - Different convergence paths

devtools::load_all()

# Parameters
d <- 5
T_per_period <- c(1000, 1000)

phi1 <- matrix(0, d, d)
diag(phi1) <- 0.3
phi1[2, 3] <- 0.4

phi2 <- matrix(0, d, d)
diag(phi2) <- 0.3
phi2[1, 2] <- 0.4

# Simulate data
set.seed(123)
sim <- var_sim_tvp(
  n = sum(T_per_period),
  sigma = diag(d) * 0.5,
  phi_list = list(phi1, phi2),
  T_per_period = T_per_period,
  type = "piecewise"
)

# Use fixed lambda for all models
fixed_lambda <- 0.1

# Extract period data
period1_data <- sim$data[[1]][1:T_per_period[1], , drop = FALSE]
period2_data <- sim$data[[1]][(T_per_period[1] + 1):sum(T_per_period), , drop = FALSE]

# Fit period-specific non-TVP models with fixed lambda
object_p1 <- constructModel(
  list(period1_data),
  tvp = FALSE,
  standardize = FALSE,
  lassotype = "standard",
  lambda1 = fixed_lambda
)
fit_p1 <- cv.multivar(object_p1)

object_p2 <- constructModel(
  list(period2_data),
  tvp = FALSE,
  standardize = FALSE,
  lassotype = "standard",
  lambda1 = fixed_lambda
)
fit_p2 <- cv.multivar(object_p2)

# Fit K=1 TVP model with common_effects=FALSE and same fixed lambda
object_tvp <- constructModel(
  sim$data,
  tvp = TRUE,
  breaks = list(sim$breaks),
  common_effects = FALSE,
  standardize = FALSE,
  lassotype = "standard",
  lambda1 = fixed_lambda
)
fit_tvp <- cv.multivar(object_tvp)

# Extract transition matrices
mat_p1 <- fit_p1$mats$total[[1]]
mat_p2 <- fit_p2$mats$total[[1]]
mat_tvp_p1 <- fit_tvp$mats$tvp[[1]][[1]]
mat_tvp_p2 <- fit_tvp$mats$tvp[[1]][[T_per_period[1]]]

# Compare
#cat("=== Period 1 Comparison (lambda =", fixed_lambda, ") ===\n")
#cat("Period-specific model:\n")
#print(round(mat_p1, 6))
#cat("\nTVP model (period 1):\n")
#print(round(mat_tvp_p1, 6))
diff_p1 <- max(abs(mat_p1 - mat_tvp_p1))
#cat("\nMax absolute difference:", diff_p1, "\n")

#cat("\n=== Period 2 Comparison (lambda =", fixed_lambda, ") ===\n")
#cat("Period-specific model:\n")
#print(round(mat_p2, 6))
#cat("\nTVP model (period 2):\n")
#print(round(mat_tvp_p2, 6))
diff_p2 <- max(abs(mat_p2 - mat_tvp_p2))
#cat("\nMax absolute difference:", diff_p2, "\n")

# Test for near-equality (within numerical tolerance)
# Use 0.01 tolerance to account for numerical differences in optimization
tol <- 0.01
p1_equal <- diff_p1 < tol
p2_equal <- diff_p2 < tol

# cat("\n=== Test Results ===\n")
# cat("Period 1 matrices nearly identical (tol =", tol, "):", p1_equal, "\n")
# cat("Period 2 matrices nearly identical (tol =", tol, "):", p2_equal, "\n")

# if (p1_equal && p2_equal) {
#   cat("\nPASS: All transition matrices are nearly identical.\n")
# } else {
#   cat("\nFAIL: Transition matrices differ more than expected.\n")
# }

expect_equal(
  p1_equal && p2_equal, TRUE
)


# Additional diagnostic: check correlation and Frobenius norm
# cat("\n=== Diagnostic Information ===\n")
# cat("Period 1 - Correlation:", cor(as.vector(mat_p1), as.vector(mat_tvp_p1)), "\n")
# cat("Period 2 - Correlation:", cor(as.vector(mat_p2), as.vector(mat_tvp_p2)), "\n")
# cat("Period 1 - Relative Frobenius:", norm(mat_p1 - mat_tvp_p1, "F") / norm(mat_p1, "F"), "\n")
# cat("Period 2 - Relative Frobenius:", norm(mat_p2 - mat_tvp_p2, "F") / norm(mat_p2, "F"), "\n")

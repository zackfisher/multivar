# Single run test with RIDGE initial estimates
# Compare Formula A vs B

library("multivar")

set.seed(123)

d <- 4
k_per_group <- 3
k <- k_per_group * 2
n <- 200

cat("=== Test with RIDGE initial estimates ===\n\n")

# Create TRUE effects with structure
true_common <- matrix(0, d, d)
true_common[1, 2] <- 0.3
true_common[2, 3] <- 0.25
true_common[3, 4] <- 0.2
diag(true_common) <- c(0.4, 0.35, 0.3, 0.25)

true_subgroup_1 <- matrix(0, d, d)
true_subgroup_1[1, 3] <- 0.15
true_subgroup_1[2, 4] <- 0.1

true_subgroup_2 <- matrix(0, d, d)
true_subgroup_2[3, 1] <- 0.12
true_subgroup_2[4, 2] <- 0.18

true_unique <- list()
set.seed(456)
for (i in 1:k) {
  u <- matrix(0, d, d)
  idx <- sample(which(true_common == 0 &
    (if(i <= k_per_group) true_subgroup_1 else true_subgroup_2) == 0), 2)
  u[idx] <- runif(2, 0.05, 0.15)
  true_unique[[i]] <- u
}

true_total <- lapply(1:k, function(i) {
  subgrp <- if (i <= k_per_group) true_subgroup_1 else true_subgroup_2
  true_common + subgrp + true_unique[[i]]
})

# Simulate data
sigma <- diag(d)
set.seed(789)
simulated_data <- lapply(1:k, function(i) {
  multivar:::var_sim(n, true_total[[i]], sigma)
})
simulated_data <- lapply(simulated_data, function(df) {
  colnames(df) <- paste0("V", 1:ncol(df))
  df
})

subgroup_membership <- c(rep(1, k_per_group), rep(2, k_per_group))

cat("True structure:\n")
cat("- Common: 7 non-zero cells (diagonal + 3 off-diagonal)\n")
cat("- Subgroup 1: 2 non-zero cells\n")
cat("- Subgroup 2: 2 non-zero cells\n")
cat("- Unique: 2 non-zero cells per subject\n\n")

# Run model with RIDGE for initial estimates
cat("Running with weightest = 'ridge'...\n")
object <- multivar::constructModel(simulated_data,
                                    subgroup_membership = subgroup_membership,
                                    weightest = "ridge")
fit <- multivar::cv.multivar(object)

initcoefs <- fit$obj@initcoefs

# Create sim_obj for eval function
sim_obj <- list(
  mat_com = true_common,
  mat_ind_unique = true_unique,
  mat_ind_final = true_total
)

cat("\n=== eval_multivar_performance results ===\n\n")
perf <- multivar::eval_multivar_performance(sim_obj, fit, averages.only = TRUE, reduced.output = TRUE)
print(perf)

cat("\n=== Initcoefs comparison ===\n")
cat("\nInitcoefs unique vs true unique (avg Frobenius):\n")
err_init_unique <- mean(sapply(1:k, function(i) {
  norm(initcoefs$unique_effects[[i]] - true_unique[[i]], "F")
}))
cat("  ", round(err_init_unique, 4), "\n")

cat("\nInitcoefs common vs true common (Frobenius):\n")
cat("  ", round(norm(initcoefs$common_effects - true_common, "F"), 4), "\n")

cat("\n=== Reconstruction check ===\n")
reconstructed <- initcoefs$common_effects + initcoefs$subgroup_effects[[1]] + initcoefs$unique_effects[[1]]
diff <- norm(initcoefs$total_effects[[1]] - reconstructed, "F")
cat("common + subgroup + unique = total? Error:", round(diff, 6), "\n")

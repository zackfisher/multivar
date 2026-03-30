# multivar

Penalized estimation of multiple-subject vector autoregressive (VAR) models. Estimates shared (common) and subject-specific (unique) network dynamics across multiple individuals using structured penalties with FISTA optimization.

## Installation

```r
# From GitHub
devtools::install_github("zackfisher/multivar")
```

## Quick start

```r
library(multivar)

# Simulate data: 2 subjects, 5 variables, 50 timepoints
sim <- multivar_sim(k = 2, d = 5, n = 50,
                    prop_fill_com = 0.1, prop_fill_ind = 0.1,
                    lb = 0.1, ub = 0.5, sigma = diag(5))

# Fit model
model <- constructModel(data = sim$data)
fit <- cv.multivar(model)

# View estimated dynamics
print_dynamics(fit)
```

## Features

- Adaptive LASSO with debiased initial estimation
- Cross-validation and eBIC model selection
- Subgroup detection
- Parallel cross-validation

## References

Fisher, Z. F., Kim, Y., Fredrickson, B. L., & Pipiras, V. (2022). Penalized estimation and forecasting of multiple subject intensive longitudinal data. *Psychometrika*, 87(4), 1377–1404.

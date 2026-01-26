## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----load-package-------------------------------------------------------------
library(multivar)

## ----fit-simple-model, results='hide', message=FALSE--------------------------
# Simulate data with K=2 subjects
set.seed(123)
sim <- multivar_sim(
  k = 2,
  d = 4,
  n = 100,
  prop_fill_com = 0.3,
  prop_fill_ind = 0.15,
  sigma = diag(4),
  lb = 0.1,
  ub = 0.9
)

# Fit the model
object <- constructModel(sim$data, intercept = FALSE)
fit <- cv.multivar(object)

## ----summary-basic------------------------------------------------------------
summary_multivar(fit)

## ----print-dynamics-fit-------------------------------------------------------
print_dynamics(fit, digits = 3)

## ----print-dynamics-custom, eval=FALSE----------------------------------------
# # More decimal places
# print_dynamics(fit, digits = 4)
# 
# # Show zeros as "0" instead of "."
# print_dynamics(fit, zero_char = "0")
# 
# # Without time subscripts
# print_dynamics(fit, time_labels = FALSE)

## ----fit-tvp-model, results='hide', message=FALSE-----------------------------
# Create time-varying dynamics
phi_period1 <- matrix(c(
  0.3,  0.4,  0.0,  0.0,
  0.0,  0.5,  0.3,  0.0,
  0.4,  0.0,  0.2,  0.0,
  0.0,  0.3,  0.0,  0.4
), nrow = 4, byrow = TRUE)

phi_period2 <- matrix(c(
  0.4,  0.3,  0.0,  0.2,
  0.0,  0.6,  0.4,  0.0,
  0.3,  0.0,  0.3,  0.0,
  0.0,  0.4,  0.0,  0.5
), nrow = 4, byrow = TRUE)

phi_period3 <- matrix(c(
  0.5,  0.2,  0.0,  0.0,
  0.0,  0.4,  0.5,  0.0,
  0.5,  0.0,  0.4,  0.2,
  0.0,  0.5,  0.0,  0.3
), nrow = 4, byrow = TRUE)

# Generate data
phi_list <- c(
  rep(list(phi_period1), 80),
  rep(list(phi_period2), 80),
  rep(list(phi_period3), 80)
)

data1 <- multivar:::var_sim_growth(
  n = 240,
  phi = phi_list,
  sigma = diag(4),
  intercept = rep(0, 4)
)

# Fit TVP model
object_tvp <- constructModel(
  list(data1),
  intercept = FALSE,
  tvp = TRUE,
  breaks = list(c(80, 160)),
  cv = "blocked",
  nfolds = 5
)

fit_tvp <- cv.multivar(object_tvp)

## ----summary-tvp--------------------------------------------------------------
summary_multivar(fit_tvp)

## ----print-dynamics-tvp-------------------------------------------------------
print_dynamics(
  fit_tvp,
  period_labels = c("Early", "Middle", "Late"),
  digits = 3
)

## ----edge-prevalence-proportion-----------------------------------------------
print_edge_prevalence(fit, type = "proportion", digits = 2)

## ----edge-prevalence-count----------------------------------------------------
print_edge_prevalence(fit, type = "count")

## ----fit-k3-model, results='hide', message=FALSE------------------------------
# Simulate K=3 model
sim3 <- multivar_sim(
  k = 3,
  d = 4,
  n = 100,
  prop_fill_com = 0.3,
  prop_fill_ind = 0.15,
  sigma = diag(4),
  lb = 0.1,
  ub = 0.9
)

object3 <- constructModel(sim3$data, intercept = FALSE)
fit3 <- cv.multivar(object3)

## ----edge-prevalence-k3-------------------------------------------------------
print_edge_prevalence(fit3, type = "proportion", digits = 2)

## ----complete-workflow, eval=FALSE--------------------------------------------
# # 1. Fit your model
# fit <- cv.multivar(object)
# 
# # 2. Get comprehensive summary
# summary_multivar(fit)
# 
# # 3. Visualize dynamics
# print_dynamics(fit)
# 
# # 4. For K>1: Check edge prevalence
# print_edge_prevalence(fit, type = "proportion")

## ----customization-examples, eval=FALSE---------------------------------------
# # Summary (no customization needed - automatic)
# summary_multivar(fit)
# 
# # Dynamics with custom formatting
# print_dynamics(
#   fit,
#   digits = 4,              # More precision
#   zero_char = "—",         # Custom zero character
#   period_labels = c(...)   # For TVP models
# )
# 
# # Edge prevalence with custom display
# print_edge_prevalence(
#   fit,
#   type = "count",          # Show counts instead of proportions
#   digits = 3,              # Decimal places for proportions
#   zero_char = "—"          # Custom zero character
# )


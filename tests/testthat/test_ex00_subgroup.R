  library(multivar)
  
  set.seed(123)
  
  #-------------------------------------------------------#
  # Parameters
  #-------------------------------------------------------#
  
  k <- 6        # Total number of subjects
  d <- 5        # Number of variables
  n <- 150      # Number of timepoints per subject
  
  # Define subgroup membership: 3 subjects in group 1, 3 in group 2
  subgroup <- c(1, 1, 1, 2, 2, 2)
  
  #-------------------------------------------------------#
  # Step 1: Simulate subgrouping data
  #-------------------------------------------------------#
  sim <- multivar_sim_subgroups(
    k = k,
    d = d,
    n = n,
    subgroup = subgroup,
    p_com = 0.15,      # Proportion of common effects (diagonal)
    p_sub = 0.10,      # Proportion of subgroup effects per group
    p_ind = 0.05,      # Proportion of unique effects per individual
    sigma = diag(d),
    lb = 0.2,          # Lower bound for coefficients
    ub = 0.6,          # Upper bound for coefficients
    max_eig_threshold = 0.95
  )
  
  #-------------------------------------------------------#
  # Step 2: Construct and fit the model
  #-------------------------------------------------------#
  object <- constructModel(
    data = sim$data,
    subgroup_membership = subgroup,   # Key argument for subgrouping
    lambda_choice = "lambda.1se",
    nfolds = 5,
    depth = 500
  )
  
  fit <- cv.multivar(object)

#saveRDS(fit, file = "/Users/zacharyfisher/Dropbox/GitHub/multivar2/tests/testthat/rds/test00_subgroup.rds")
  
  mats <- readRDS("rds/test00_subgroup.rds")$mats

  #-------------------------------------------------------# 
  context("test00: common effects correct")
  #-------------------------------------------------------# 
  
  expect_equal_to_reference(
    fit$mats$common, mats$common
  )
  
  #-------------------------------------------------------# 
  context("test00: subgroup effects correct")
  #-------------------------------------------------------# 
  
  #mats <- readRDS("rds/test00_common_effects.rds")
  
  expect_equal_to_reference(
    fit$mats$subgrp, mats$subgrp
  )
  
  #-------------------------------------------------------# 
  context("test00: unique effects correct")
  #-------------------------------------------------------# 
  
  expect_equal_to_reference(
    fit$mats$unique, mats$unique
  )
  
  #-------------------------------------------------------# 
  context("test00: total effects correct")
  #-------------------------------------------------------# 
  
  expect_equal_to_reference(
    fit$mats$total, mats$total
  )
library("multivar")

  set.seed(123)
  
  d <- 5
  k <- 3
  n <- 100
  
  sim  <- multivar::multivar_sim(
    k = k,   # number of individuals
    d = d,   # number of variables
    n = n,   # number of timepoints
    prop_fill_com = .2,
    prop_fill_ind = .1,
    unique_overlap = TRUE,
    sigma = diag(d),
    lb = .1,
    ub =.9
  )
  
  object <- multivar::constructModel(sim$data, lambda_choice = "lambda.1se", nfolds = 10,depth=1000)
  fit <- multivar::cv.multivar(object)
  
  #saveRDS(fit$mats$common, file = "tests/testthat/rds/test00_common_effects.rds")
  #saveRDS(fit$mats$unique, file = "tests/testthat/rds/test00_unique_effects.rds")
  #saveRDS(fit$mats$total,  file = "tests/testthat/rds/test00_total_effects.rds")
  #old_com <- readRDS("tests/testthat/rds/test00_common_effects.rds")
  #fit$mats$common
  #-------------------------------------------------------# 
  context("test00: common effects correct")
  #-------------------------------------------------------# 
  
  expect_equal_to_reference(
    fit$mats$common, "rds/test00_common_effects.rds"
  )
  
  #-------------------------------------------------------# 
  context("test00: unique effects correct")
  #-------------------------------------------------------# 
  
  expect_equal_to_reference(
    fit$mats$unique, "rds/test00_unique_effects.rds"
  )
  
  #old_unique <- readRDS("tests/testthat/rds/test00_unique_effects.rds")
  # old_unique[[1]]
  # fit$mats$unique[[1]]
  # sim$mat_ind_unique[[1]]
  # 
  # old_unique[[4]]
  # fit$mats$unique[[4]]
  # sim$mat_ind_unique[[4]]
  
  #-------------------------------------------------------# 
  context("test00: total effects correct")
  #-------------------------------------------------------# 
  
  expect_equal_to_reference(
    fit$mats$total, "rds/test00_total_effects.rds"
  )
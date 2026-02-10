get_subgroups <- function(data, nlambda1, pendiag){

  # Construct and fit no subgrouping model
  no_sub    <- constructModel(data = data, nlambda1 = nlambda1, pendiag = pendiag)
  fit_nosub <- cv.multivar(no_sub) 
  
  # Get subgroup membership
  k       <- no_sub@k
  mats    <- fit_nosub$mats$unique
  sim_vec <- c()
  
  for (i in 1:(length(mats)-1)){
    for (j in (i+1):length(mats)){
      sim_count <- sum( (mats[[i]] > 0 & mats[[j]] > 0) | (mats[[i]] < 0 & mats[[j]] < 0) )
      sim_vec <- c(sim_vec, sim_count)
    }
  }
  
  sim_mat <- matrix(0, nrow = k, ncol = k)
  sim_mat[lower.tri(sim_mat)] <- sim_vec
  sim_mat <- sim_mat + t(sim_mat)
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for get_subgroups(). Please install it.", call. = FALSE)
  }
  
  g <- igraph::graph_from_adjacency_matrix(
    sim_mat,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  
  weights <- igraph::E(g)$weight
  res     <- igraph::cluster_walktrap(g, weights = weights, steps = 4)
  sub_mem <- as.numeric(igraph::membership(res))
  
  return(sub_mem)
}

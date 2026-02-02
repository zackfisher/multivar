#' Simulate multivar data with time-varying dynamics according to a latent growth curve.
#'
#' @param k Integer. The number of individuals (or datasets) to be generated.
#' @param d Integer. The number of variables per dataset. For now this will be constant across individuals. 
#' @param n Integer. The time series length. 
#' @param prop_fill_com Numeric. The proportion of nonzero paths in the common transition matrix.
#' @param prop_fill_ind Numeric. The proportion of nonzero unique (not in the common transition matrix or transition matrix of other individuals) paths in each individual transition matrix.
#' @param prop_fill_com_growth Numeric. The proportion of nonzero time-varying paths in the common transition matrix.
#' @param prop_fill_ind_growth Numeric. The proportion of nonzero time-varying unique (not in the common transition matrix or transition matrix of other individuals) paths in each individual transition matrix.
#' @param offset Numeric. An offset dictating when dynamics should transition from time-invariant to time-varying.
#' @param ub Numeric. The lower bound for individual elements of the transition matrices.
#' @param lb Numeric. The upper bound for individual elements of the transition matrices.
#' @param sigma Matrix. The (population) innovation covariance matrix.
#' @param diag Logical. Default is FALSE. Should diagonal elements be filled first for common elements.
#' @param intercept List. Default is NULL. A list of length K containing numeric vectors of length d representing the intercept values. If NULL, intercepts are set to 0.
#' @keywords var multivar simulate
#' @export
multivar_sim_breaks <- function(
  k,
  d,
  n,
  prop_fill_com,
  prop_fill_ind,
  prop_fill_com_growth,
  prop_fill_ind_growth,
  offset = 0,
  lb,
  ub,
  lb_int,
  ub_int,
  lb_slope,
  ub_slope,
  sigma,
  diag = FALSE,
  iters = 1,
  intercept = NULL){

  # k = 30
  # d = 10
  # n = 50
  # prop_fill_com = 0.1
  # prop_fill_ind = 0.03
  # prop_fill_com_growth  = 0.03
  # prop_fill_ind_growth  = 0.03
  # offset = 0
  # lb = 0.1
  # ub = 0.9
  # lb_int = 0
  # ub_int = 0
  # lb_slope = .8/n
  # ub_slope = .8/n
  # sigma = diag(d)
  # unique_overlap = TRUE
  # mat_common = NULL
  # mat_unique = NULL
  # mat_total = NULL
  # diag = TRUE
  
  d.fun                 <- function(n) runif(n,lb,ub)
  densities             <- c(prop_fill_com, rep(prop_fill_ind,k)) 
  
  n_elem_prop_fill_com  <- round(d^2*prop_fill_com) 
  n_elem_prop_fill_ind  <- round(d^2*prop_fill_ind) 
  n_elem_total_all_subj <- n_elem_prop_fill_ind*k
  
  n_elem_prop_fill_com_growth  <- round(d^2*prop_fill_com_growth) 
  n_elem_prop_fill_ind_growth  <- round(d^2*prop_fill_ind_growth) 
  n_elem_total_all_subj_growth <- n_elem_prop_fill_ind_growth*k
  
  n_elem_total_all <- n_elem_total_all_subj + 
                      n_elem_prop_fill_com + 
                      n_elem_total_all_subj_growth + 
                      n_elem_prop_fill_com_growth
  
  max_eigs <- big_values <- TRUE
  
  # Set default intercepts if not provided
  if(is.null(intercept)){
    intercept <- replicate(k, rep(0, d), simplify = FALSE)
  }

  # the order of this matters
  props <- c(
    prop_fill_com,
    prop_fill_com_growth,
    prop_fill_ind,
    prop_fill_ind_growth
  )

  while(max_eigs | big_values){
  
    M <- array(0, dim = list(d^2, k, n))
    
    path_df <- data.frame()
    
    bad_paths <- replicate(k,c(),simplify=FALSE)
      
    A  <- matrix(c(1:d^2),d,d)
    
    # if diag_first we want the diagonal paths to go first
    if(diag){
      vec_A             <- c(do.call(rbind,lapply(1:nrow(A),function(r){setdiff(sample(A[r,]),diag(A)[r])})))
      vec_A_diag        <- diag(matrix(c(1:d^2),d,d)) # indices for diagonal elements of A
      vec_A_nondiag     <- setdiff(vec_A,vec_A_diag)  # indices for non-diagonal elements of A
      permuted_elements <- c(vec_A_diag,vec_A_nondiag)
    } else {
      permuted_elements <- c(do.call(rbind,lapply(1:nrow(A),function(r){sample(A[r,])})))
    }
    #------------------------------------------------------------------------#
    # common effects (time-invariant)
    #------------------------------------------------------------------------#
    if(props[1] != 0){ # j <- 1
      path_idx <- setdiff(permuted_elements,bad_paths)[1:round(props[1]*d^2)]
      subj_idx <- c(1:k)
      m_idx <- expand.grid(path_idx,subj_idx)
      m_idx <- cbind(m_idx, runif(nrow(m_idx), lb,ub)) 

      for(i in 1:dim(M)[3]){
        M0 <- matrix(0, d^2, k)
        M0[as.matrix(m_idx[1:2])] <- m_idx[,3]
        M[,,i] <- M0
      }
      bad_paths <- lapply(bad_paths, function(x){path_idx})
      
      m_1 <- cbind("com",m_idx,matrix(m_idx[,3],nrow(m_idx),n))
      colnames(m_1) <- c("type","elem","subj","int",paste0("t",1:n))
    }
      
    #------------------------------------------------------------------------#
    # common effects (time-varying)
    #------------------------------------------------------------------------#
    if(props[2] != 0){ # 
      
      path_idx <- setdiff(permuted_elements,bad_paths[[1]])[1:round(props[2]*d^2)]
      subj_idx <- c(1:k)
      m_idx <- expand.grid(path_idx,subj_idx)
  
      # generate intercept and slope 
      G <- matrix(0, k*length(path_idx), n + 1)
      
      if(offset != 0){
        
        G[,c(1:(offset+1))] <- runif(1, lb_int,ub_int)
        
        for(i in 1:nrow(G)){
          
          G[i,c((offset+2):ncol(G))] <- G[i,1] + runif(1,lb_slope,ub_slope)*c(1:(n-offset))
          
        }
        
      } else {
        
        G[,1] <- runif(1, lb_int,ub_int)
        
        for(i in 1:nrow(G)){
          
          G[i,2:ncol(G)] <- G[i,1] + runif(1,lb_slope,ub_slope)*c(1:n)
          
        }
        
      }
      
      m_idx <- cbind(m_idx, G) 
      colnames(m_idx) <- c("elem","subj","int",paste0("t",1:n))
      
      for(i in 1:dim(M)[3]){
        M0 <- M[,,i]
        M0[as.matrix(m_idx[1:2])] <- m_idx[,3+i]
        M[,,i] <- M0
      }
      bad_paths <- lapply(bad_paths,function(x) {c(x,path_idx)})
      
      m_idx <- cbind("com_tvp",m_idx)
      colnames(m_idx)[1] <- "type"
      
      m_2 <- rbind(m_1, m_idx)
      
    } else {
      
      m_2 <- m_1
      
    }
    
    #------------------------------------------------------------------------#
    # unique effects (time-invariant)
    #------------------------------------------------------------------------#
    if(props[3] != 0){

      path_idx <- list()
      for(i in 1:k){
        poss_paths <- c(do.call(rbind,lapply(1:nrow(A),function(r){setdiff(sample(A[r,]),diag(A)[r])})))
        path_idx[[i]]  <- setdiff(poss_paths,bad_paths[[i]])[1:round(props[3]*d^2)]
        bad_paths[[i]] <- c(bad_paths[[i]],path_idx[[i]])
        #path_idx[[i]]  <- sample(setdiff(permuted_elements,bad_paths[[i]]),round(props[3]*d^2))
        #bad_paths[[i]] <- c(bad_paths[[i]],path_idx[[i]])
      }

      m_idx <- do.call(rbind,lapply(seq_along(path_idx),function(kk){cbind(path_idx[[kk]],kk)}))
      m_idx <- cbind(m_idx, runif(nrow(m_idx), lb,ub))

      for(i in 1:dim(M)[3]){
        M0 <- M[,,i]
        M0[as.matrix(m_idx[,1:2])] <- m_idx[,3]
        M[,,i] <- M0
      }
      m_idx <- as.matrix(m_idx)
      m_3 <- cbind("ind",cbind(m_idx,matrix(m_idx[,3],nrow(m_idx),n)))
      colnames(m_3) <- c("type","elem","subj","int",paste0("t",1:n))
      m_3 <- rbind(m_2,m_3)
    } else {
      m_3 <- m_2
    }
    
    #------------------------------------------------------------------------#
    # unique effects (time-varying)
    #------------------------------------------------------------------------#
    if(props[4] != 0){
      
      path_idx <- list()
      for(i in 1:k){
        poss_paths <- c(do.call(rbind,lapply(1:nrow(A),function(r){setdiff(sample(A[r,]),diag(A)[r])})))
        path_idx[[i]]  <- setdiff(poss_paths,bad_paths[[i]])[1:round(props[4]*d^2)]
        bad_paths[[i]] <- c(bad_paths[[i]],path_idx[[i]])
        #path_idx[[i]]  <- sample(setdiff(permuted_elements,bad_paths[[i]]),round(props[4]*d^2))
        #bad_paths[[i]] <- c(bad_paths[[i]],path_idx[[i]])
      }
      
      m_idx <- do.call(rbind,lapply(seq_along(path_idx),function(kk){cbind(path_idx[[kk]],kk)}))
      
      # generate intercept and slope 
      G <- matrix(0, k*length(path_idx[[1]]), n + 1)
      
      if(offset != 0){
        
        G[,c(1:(offset+1))] <- runif(1, lb_int,ub_int)
        
        for(i in 1:nrow(G)){
          
          G[i,c((offset+2):ncol(G))] <- G[i,1] + runif(1,lb_slope,ub_slope)*c(1:(n-offset))
          
        }
        
      } else {
        
        G[,1] <- runif(1, lb_int,ub_int)
        
        for(i in 1:nrow(G)){
          
          G[i,2:ncol(G)] <- G[i,1] + runif(1,lb_slope,ub_slope)*c(1:n)
          
        }
        
      }
      
      m_idx <- cbind(m_idx, G) 
      colnames(m_idx) <- c("elem","subj","int",paste0("t",1:n))
      
      for(i in 1:dim(M)[3]){
        M0 <- M[,,i]
        M0[as.matrix(m_idx[,1:2])] <- m_idx[,3+i]
        M[,,i] <- M0
      }
      
      m_idx <- cbind("ind_tvp",m_idx)
      colnames(m_idx)[1] <- "type"

      m_4 <- rbind(m_3, m_idx)

    } else {

      m_4 <- m_3

    }


    #------------------------------------------------------------------------#
    # assemble transition mats
    #------------------------------------------------------------------------#
    phi <- replicate(n, matrix(0,d,d), simplify = FALSE)
    phi <- replicate(k, phi, simplify = FALSE)
    for(kk in 1:k){
      for(nn in 1:n){
        phi[[kk]][[nn]] <- matrix(M[,kk,nn],d,d)
      }
    }
        
    max_eigs <- any(unlist(
      lapply(phi, function(mats){
        unlist(lapply(mats, function (x){
          max(abs(eigen(x)$values))
        }))
      })
    ) > .99)
    
    # while(max_eigs){
    #   
    #   for(kk in 1:k){
    #     for(nn in 1:n){
    #       phi[[kk]][[nn]] <- phi[[kk]][[nn]]*.99
    #     }
    #   }
    #   
    #   max_eigs <- any(unlist(
    #     lapply(phi, function(mats){
    #       unlist(lapply(mats, function (x){
    #         max(abs(eigen(x)$values))
    #       }))
    #     })
    #   ) > .99)
    #   
    # }
    
    # transition_mats <- lapply(transition_mats, function(x){
    #   max_eig <- max(abs(eigen(x)$values))
    #   if(max_eig > .99){
    #     x <- x*(.9/max_eig)
    #   }
    #   x
    # })
    
    cnter <- 1
    
    if(!max_eigs){
    
      while (big_values & cnter < 100){
      
        data <- list()

        for(jj in 1:iters){
          #print(paste0("iter: ", jj))
          data_jj <- lapply(seq_along(phi), function(i) {var_sim_growth(n, phi[[i]], sigma, intercept = intercept[[i]])})
          data_jj <- lapply(data_jj, function(df){colnames(df) <- paste0("V",1:ncol(df)); df})
          data[[jj]] <- data_jj
        }
        
        big_values <- any(unlist(lapply(data, function(df){
          any(unlist(lapply(seq_along(df),function(kk){any(abs(df[[kk]])>15)})))
        })))
        
        cnter <- cnter + 1
        #print(cnter)
      
      }
    
    #if(big_values){ print(paste0("big_values"))}
    #if(max_eigs){ print(paste0("max_eigs"))}
    }
    
  }
  
  phi_types <- lapply(1:k, function(kk){
    tc <- m_4[m_4$subj==kk,c("type","elem")]
    tb <- rep(0,d^2)
    tb[as.numeric(tc$elem)] <- tc$type
    ta <- matrix(tb, d, d)
    ta
  })
    
  mat_list <- list(
    mat_com        = NULL,
    mat_ind_unique = NULL,
    mat_ind_final  = lapply(1:k, function(i){phi[[i]]}),
    data = data,
    path_idx = m_4,
    offset = offset,
    phi = phi,
    phi_types = phi_types,
    intercept = intercept
  )

  return(mat_list)

}


fit_canonical_var <- function(A, p = 1, type = "none"){

    fit <- vars::VAR(A, p = p, type = type)

    d <- ncol(A)

    # Extract transition matrix (Phi)
    # If type="const", coefficients include intercept as last element
    transition_mat <- as.matrix(do.call("rbind",lapply(1:d, function(x) {
      coefs <- fit$varresult[[x]]$coefficients
      if (type == "const") {
        coefs[1:d]  # Exclude intercept (last element)
      } else {
        coefs
      }
    })))
    colnames(transition_mat) <-  rownames(transition_mat) <- colnames(A)

    # Extract intercepts if type="const"
    intercepts <- NULL
    if (type == "const") {
      intercepts <- sapply(1:d, function(x) {
        coefs <- fit$varresult[[x]]$coefficients
        coefs[d + 1]  # Intercept is after the d lag coefficients
      })
      names(intercepts) <- colnames(A)
    }

    transition_mat_pval <- as.matrix(do.call("rbind",lapply(1:d, function(x) {
       pvals <- coef(fit)[[x]][,"Pr(>|t|)"]
       if (type == "const") {
         pvals[1:d]  # Exclude intercept row
       } else {
         pvals
       }
    })))
    colnames(transition_mat_pval) <-  rownames(transition_mat_pval) <- colnames(A)

    transition_mat_thresh <- as.matrix(do.call("rbind",lapply(1:d, function(x) {
      pvals <- coef(fit)[[x]][,"Pr(>|t|)"]
      if (type == "const") {
        pvals[1:d] < 0.05  # Exclude intercept row
      } else {
        pvals < 0.05
      }
    })))
    colnames(transition_mat_thresh) <-  rownames(transition_mat_thresh) <- colnames(A)

    transition_mat_sigonly <- matrix(0, nrow(transition_mat), ncol(transition_mat))
    transition_mat_sigonly[transition_mat_thresh] <- transition_mat[transition_mat_thresh]
    colnames(transition_mat_sigonly) <-  rownames(transition_mat_sigonly) <- colnames(A)

    return(list(
      transition_mat = transition_mat,
      transition_mat_sigonly = transition_mat_sigonly,
      intercepts = intercepts
    ))


}
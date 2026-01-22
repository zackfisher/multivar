#' @export
setup_data <- function (data, standardize, lag, horizon, intercept) {
  
   if (is.null(data)){
    stop(paste0(
      "multivar ERROR: neither a data directory nor a data list is specified. ",
      "Please either specify a directory of data or a list of individual data files."
    )
    )
  }
  
  #-------------------------------------------------------------#
  # If the data is already in list form.
  #-------------------------------------------------------------#
  if (!is.list(data)){
      stop(paste0(
        "multivar ERROR: data must be supplied as a list of matrices."
      ))
  } else {
      ts_list  <- list()
  
      # if the user-supplied list does not have names, add names
      if(is.null(names(data))){
  
        names(data) <- paste0("dataset", 1:length(data))
  
      }
  
      ts_list  <- data
  }

  
  #-------------------------------------------------------------#
  # Ensure all datafiles share the same column order.
  #
  #   # Now we also make all datasets have the same column order.
  #-------------------------------------------------------------# 
  varnames       <- colnames(ts_list[[1]])
  n_orig_vars    <- ncol(ts_list[[1]])
  
  if(is.null(varnames)){
    varnames <- c(paste0("V", seq(1,n_orig_vars)))
    ts_list <- lapply(ts_list, function(x) { colnames(x) <- varnames;  x })
  } 
  
  
  #-------------------------------------------------------------#
  # Final data manipulation.
  #-------------------------------------------------------------#
  if (ncol(ts_list[[1]]) == 1) {
    stop(paste0("multivar ERROR: only one column of data read in. ",
                "Check if sep argument properly specified."))
  }
  
  

  if(standardize & intercept == FALSE){
    ts_list <- lapply(ts_list, function(df) {scale(df)})
  } else if(standardize & intercept == TRUE){
    ts_list <- lapply(ts_list, function(df) {scale(df, center = FALSE)})
  }

  
  #-------------------------------------------------------------#
  # Final data checks
  #-------------------------------------------------------------#
  
  n_subjects   <- length(ts_list)
  cols         <- numeric()
  missingData  <- numeric()
  constantCols <- logical()
  numericCols  <- logical()
  
  # check for obvious errors in data
  for (k in 1:length(ts_list)){
    data.file <- ts_list[[k]]
    cols[k]   <- ncol(data.file)
    missingData[k]  <- any(is.na(data.file))
    constantCols[k] <- any(apply(data.file, 2, sd, na.rm = TRUE) == 0)
    numericCols[k]  <- any(apply(data.file, 2, is.numeric) == FALSE)
  }
  
  # will need to be updated for lags > 1
  ts_list <- lapply(ts_list, function(df){
    if(horizon > 0){
      H  <- df[((nrow(df)-horizon+1):nrow(df)),, drop = FALSE]
      df <- df[-((nrow(df)-horizon+1):nrow(df)),, drop = FALSE]
    } else {
      H  <- NA
    }
    A  <- Matrix(df[1:(nrow(df)-1), ], sparse = FALSE)
    b  <- Matrix(df[2:(nrow(df)  ), ], sparse = FALSE)
    colnames(A) <- colnames(b) <- colnames(df)
    list(b = b, A = A, H = H)
  })
  
  if (intercept){
    ts_list <- lapply(ts_list, function(lst){
      A_int <- cbind(1, lst$A)
      colnames(A_int)[1] <- "Intercept"
      lst$A <- A_int
      lst
    })
  }
  
  return(ts_list)
  
}
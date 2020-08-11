#' Check multivar object.
check.multivar <- function(object){

    errors <- character()

    # if(any(is.na(Y))){
    #   msg <- c("Remove NA values before running ConstructModel")
    #   errors <- c(errors,msg)
    # }      
    
    if(length(errors)==0) TRUE else errors
    
}

#' multivar Object Class
#'
#' An object class to be used with cv.multivar
#' 
#' @slot n Numeric Vector. Vector containing the number of variables for each dataset.
#' @slot d Numeric Vector. Vector containing the number of timepoints for each dataset.
#' @slot Ak List. A list (length = k) of lagged (T-lag-horizon) by d multivariate time series.
#' @slot bk List. A list (length = k) of (T-lag-horizon) by d multivariate time series.
#' @slot bk List. A list (length = k) of (horizon) by d multivariate time series.
#' @slot A  Matrix. A matrix containing the lagged ((T-lag-horizon)k) by (d+dk) multivariate time series.
#' @slot b  Matrix. A matrix containing the non-lagged ((T-lag-horizon)k) by (d) multivariate time series.
#' @slot H  Matrix. A matrix containing the non-lagged (horizon k) by d multivariate time series.
#' @slot lag Numeric. The VAR order. Currently only lag 1 is supported.
#' @slot horizon Numeric. Forecast horizon.
#' @slot t1 Numeric vector. Index of time series in which to start cross validation for individual k. 
#' @slot t2 Numeric vector. Index of time series in which to end cross validation for individual k.
#' @slot tol Numeric. Convergence tolerance.
#' @slot window Numeric. Size of rolling window.
#' @details To construct an object of class multivar, use the function \code{\link{constructModel}}
#' @seealso \code{\link{constructModel}}
#' @export
setClass(
    Class="multivar",
    representation(
        n  = "numeric",
        d  = "numeric",
        Ak ="list",
        bk = "list",
        Hk = "list",
        A  = "matrix",
        b  = "matrix",
        H  = "matrix",
        lag="numeric",
        horizon="numeric",
        t1k = "numeric",
        t2k = "numeric",
        lambdas="numeric",
        tol="numeric",
        window="numeric"
        ),validity=check.multivar
    )





#' Construct an object of class multivar
#' 
#' @param data List. A list (length = k) of T by d multivariate time series
#' @param lag Numeric. The VAR order. Default is 1.
#' @param struct Character. The choice of penalty structure can be "Sparse","Lowrank",or "SparseLowrank".
#' @param horizon Numeric. Desired forecast horizon. Default is 1.
#' @param t1 Numeric. Index of time series in which to start cross validation. If NULL, default is floor(nrow(n)/3) where nk is the time series length for individual k.
#' @param t2 Numeric. Index of times series in which to end cross validation. If NULL, default is floor(2*nrow(n)/3) where nk is the time series length for individual k.
#' @param lambdas List. User-supplied penalty parameters.
#' @param tol Numeric. Optimization tolerance (default 1e-4)
#' @param window Numeric. Size of rolling window.  If set to 0 an expanding window will be used. 
#' @param standardize Logical. Default is true. Whether to standardize the individual data.
#' 
#' @details The choices for "struct" are as follows
#' \itemize{
#' \item{  "Sparse" (Sparse)}
#' \item{  "Lowrank" (Low-Rank)} 
#' \item{  "SparseLowrank" (Sparse and Low-Rank)} 
#' }
#'
#' 
#' @references 
#' @examples
#' # multivar Example
#' @export
constructModel <- function( data = NULL,
                            lag = 1,
                            struct = "Sparse",
                            horizon=1,
                            t1 = NULL, 
                            t2 = NULL, 
                            lambdas=NULL,
                            tol=1e-4,
                            window = 1,
                            standardize = T,
                            sep = NULL,
                            header = NULL ){
  
  if( lag != 1 ){
    stop("multivar ERROR: Currently only lag of order 1 is supported.")
  }
  
  structures <- c("Sparse","Lowrank","SparseLowrank")
  if(!struct %in% structures){
    stop(cat("multivar ERROR: Struct argument must be one of the following:", structures))
  }
  
  if(!is.null(horizon) & horizon<1 ){
    stop("Forecast Horizon must be at least 1")
  }
  
  if( tol<0 | tol>1e-1 ){
    stop("Tolerance must be positive")
  }
  
  dat <- setup_data(data, sep, header, standardize, lag, horizon) 
  
  getj  <- function(mat){dp = diff(mat@p);rep(seq_along(dp), dp) - 1}
  k     <- length(dat)
  p     <- ncol(dat[[1]]$A)
  ns    <- sapply(dat, function(item){nrow(item$A)})
  cns   <- endrows <- cumsum(ns)
  sr    <- c(1,cns[-length(cns)] + 1) - 1
  nz    <- tail(cns,1)
  is    <- js <- xs <- NULL
  	
  for(ii in 1:k){
  	is <- c(is, sr[ii] + dat[[ii]]$A@i)
  	js <- c(js, getj(dat[[ii]]$A))
  	xs <- c(xs, dat[[ii]]$A@x)
  	is <- c(is, sr[ii] + dat[[ii]]$A@i)
  	js <- c(js, getj(dat[[ii]]$A) + p*ii)
  	xs <- c(xs, dat[[ii]]$A@x)
  }
  
  Ak <- lapply(dat, "[[", "A")
  bk <- lapply(dat, "[[", "b")
  Hk <- lapply(dat, "[[", "H")
  A  <- sparseMatrix(i = is, j = js, x = xs, index1=FALSE, dims = c(nz,p*(k+1)))
  b  <- as.matrix(do.call(rbind, bk))
  H  <- as.matrix(do.call(rbind, Hk))
  
  # what indices do we need for forecasting
  t1k <- unlist(lapply(dat, function(x){floor(nrow(x$b)/3)}))
  t2k <- unlist(lapply(dat, function(x){floor(2*nrow(x$b)/3)}))
  ntk <- unlist(lapply(dat, function(x){nrow(x$b)})) # number cols 
  nvk <- unlist(lapply(dat, function(x){ncol(x$b)})) # number tpts
  # tks <- c(1, cumsum(ntk[-length(ntk)])+1)
  # tke <- cumsum(ntk)
  # t1s <- c(1, cumsum(t1k[-length(t1k)])+1)
  # t1e <- cumsum(t1k)
  # t2s <- t1e+1
  # t1e <- cumsum(t1k)
  
  obj <- new("multivar",
    n  = ntk,
    d  = ndk,
    Ak = Ak,
    bk = bk,
    Hk = Hk,
    A  = A,
    b  = b,
    H  = H,
    horizon = horizon,
    t1=t1k,
    t2=t2k,
    lambdas=lambdas,
    tol=tol,
    window=window
  )

  return(obj)

}




# show-default method to show an object when its name is printed in the console.
#' Default show method for an object of class multivar
#'
#' @param object \code{multivar} object created from \code{ConstructModel}
#' @return Displays the following information about the multivar object:
#' \itemize{
#' \item{To do.}
#' }
#' @seealso \code{\link{constructModel}} 
#' @name show.multivar
#' @aliases show,multivar-method
#' @docType methods
#' @rdname show-methods
#' @export
setMethod("show","multivar",function(object){
  cat("multivar model: \n")
})


#' Cross Validation for multivar
#' 
#' @usage cv.multivar(object)
#' @param object multivar object built using \code{ConstructModel}.
#' @details The main function of the multivar package. Performs cross validation to select penalty parameters over a training sample and evaluates them over a test set. Creates an object of class \code{multivar.results}
#' @return An object of class \code{multivar.results}.
#' @seealso \code{\link{constructModel}}, \code{\link{multivar.results}},\code{\link{multivar.est}} 
#' @name cv.multivar
#' @aliases cv.multivar,multivar-method
#' @docType methods
#' @rdname cv.multivar-methods
#' @examples
#'  # Add example
#' @export
setGeneric(name="cv.multivar",def=function(object){standardGeneric("cv.multivar")})



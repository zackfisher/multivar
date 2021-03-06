% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/var_forecast.R
\name{var_forecast}
\alias{var_forecast}
\title{Estimate h-step ahead forecasts based on the recovered transition matrix.}
\usage{
var_forecast(yf, h, A)
}
\arguments{
\item{yf}{A d x T data matrix where d is the number of observed variables and T is the number of timepoints.}

\item{h}{An integer indicating the forecast horizon.}

\item{A}{A d x d transition matrix.}
}
\description{
Estimate h-step ahead forecasts based on the recovered transition matrix.
}
\examples{

theta    <- diag(c(.7,.8,.9,.6,.7,.9))
data     <- var_sim(100, theta, diag(.1,6))
datalag  <- embed(data, 2)
b        <- datalag[,1:6]
A        <- datalag[,7:12]
A_est    <- fista_sparse(A, b, 1, theta, niter = 10, backtrack = TRUE)$out.x
var_forecast(t(b), 2, A_est)

}
\keyword{forecast}
\keyword{var}

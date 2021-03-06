% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multivarObjectClass.R
\docType{class}
\name{multivar-class}
\alias{multivar-class}
\title{multivar Object Class}
\description{
An object class to be used with cv.multivar
}
\details{
To construct an object of class multivar, use the function \code{\link{constructModel}}
}
\section{Slots}{

\describe{
\item{\code{n}}{Numeric Vector. Vector containing the number of variables for each dataset.}

\item{\code{d}}{Numeric Vector. Vector containing the number of timepoints for each dataset.}

\item{\code{Ak}}{List. A list (length = k) of lagged (T-lag-horizon) by d multivariate time series.}

\item{\code{bk}}{List. A list (length = k) of (T-lag-horizon) by d multivariate time series.}

\item{\code{bk}}{List. A list (length = k) of (horizon) by d multivariate time series.}

\item{\code{A}}{Matrix. A matrix containing the lagged ((T-lag-horizon)k) by (d+dk) multivariate time series.}

\item{\code{b}}{Matrix. A matrix containing the non-lagged ((T-lag-horizon)k) by (d) multivariate time series.}

\item{\code{H}}{Matrix. A matrix containing the non-lagged (horizon k) by d multivariate time series.}

\item{\code{lag}}{Numeric. The VAR order. Currently only lag 1 is supported.}

\item{\code{horizon}}{Numeric. Forecast horizon.}

\item{\code{t1}}{Numeric vector. Index of time series in which to start cross validation for individual k.}

\item{\code{t2}}{Numeric vector. Index of time series in which to end cross validation for individual k.}

\item{\code{tol}}{Numeric. Convergence tolerance.}

\item{\code{window}}{Numeric. Size of rolling window.}
}}

\seealso{
\code{\link{constructModel}}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{confidence_band}
\alias{confidence_band}
\title{Confidence Bands of Smooth Terms}
\usage{
confidence_band(theta, V, X, q = 0.05, R = 1000)
}
\arguments{
\item{theta}{The estimated coefficients which corresponds to the smooth
term.}

\item{V}{The covariance matrix of the estimated coefficients \code{theta}.}

\item{X}{The design matrix of the model which corresponds to the smooth
term.}

\item{q}{The quantile. Default to \code{q = 0.05} which corresponds to
95\% confidence bands.}

\item{R}{The number of replications in the simulation process.}
}
\value{
A list of two vectors which contain the lower and the upper limits
of the confidence band.
}
\description{
\code{confidence_band} computes the lower and upper limits of smooth
terms fitted with \code{intensity_pspline}.
}
\keyword{internal}

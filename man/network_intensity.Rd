% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network.R
\name{network_intensity}
\alias{network_intensity}
\title{Fitted Intensity on a Geometric Network}
\usage{
network_intensity(z, m, fit1, fit2 = NULL, scale = NULL)
}
\arguments{
\item{z}{The shortest path distance from the beginning of the network
segment.}

\item{m}{The network segment index.}

\item{fit1}{A fitted geometric network.}

\item{fit2}{A second fitted geometric network. If specified, the function
returns the squared difference of the intensity fits at the specified point
of the network.}

\item{scale}{A numeric vector of length two which determines the scaling
of the two intensity functions.}
}
\value{
A numeric vector of length one, indicating the intensity (or the
squared difference of two intensities) at the specified point.
}
\description{
Fitted Intensity on a Geometric Network
}
\keyword{internal}

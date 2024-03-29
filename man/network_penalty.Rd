% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrices.R
\name{network_penalty}
\alias{network_penalty}
\title{Penalty Matrix of a Geometric Network}
\usage{
network_penalty(G, knots, r)
}
\arguments{
\item{G}{A geometric network (object of class \code{gnpp}).}

\item{knots}{A list which contains the knots on which the
B-splines are defined}

\item{r}{The order of the penalty, default to first-order penalty (\code{r = 1}.}
}
\value{
A sparse and square penalty matrix.
}
\description{
\code{network_penalty} constructs the penalty matrix which relates to the
B-Splines created by \code{\link[geonet]{bspline_design}}.
}
\author{
Marc Schneble \email{marc.schneble@stat.uni-muenchen.de}
}
\keyword{internal}

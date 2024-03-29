% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.R
\name{as_gn}
\alias{as_gn}
\alias{as_gn.linnet}
\alias{as_gn.gnpp}
\alias{as_gn.gnppfit}
\alias{as_gn.lpp}
\alias{as_gn.lppfit}
\title{Transmute to Geometric Network}
\usage{
as_gn(x, ...)

\method{as_gn}{linnet}(x, ..., spatstat = FALSE)

\method{as_gn}{gnpp}(x, ...)

\method{as_gn}{gnppfit}(x, ...)

\method{as_gn}{lpp}(x, ..., spatstat = FALSE)

\method{as_gn}{lppfit}(x, ...)
}
\arguments{
\item{x}{An object that could reasonably be transmuted to an
object of class \code{gn}.}

\item{...}{Other arguments.}

\item{spatstat}{Set to \code{TRUE} if retransformation to an object of any
\code{spatstat} class is desired. Some elements of these objects (such as
the window) are otherwise discarded when being transmuted to an object of
class \code{gn}.}
}
\value{
An object of class \code{gn}.
}
\description{
\code{as_gn} transmutes an existing object into a geometric network, an object
of class \code{gn}.
}
\examples{
library(spatstat.data)
G <- as_gn(simplenet)
summary(G)
}
\author{
Marc Schneble \email{marc.schneble@stat.uni-muenchen.de}
}

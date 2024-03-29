% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cu_crw_predict.R
\name{cu_crw_predict}
\alias{cu_crw_predict}
\title{Batch CRW Prediction for Multiple Animals}
\usage{
cu_crw_predict(
  fit,
  predTime = NULL,
  barrier = NULL,
  vis_graph = NULL,
  as_sf = TRUE,
  ...
)
}
\arguments{
\item{fit}{A CRW fit object}

\item{predTime}{A character string describing the desired frequency of prediction,
e.g., \code{predTime="1 hour"} or \code{predTime="15 min"}.}

\item{barrier}{An \code{sf} polygon object representing areas where the animal cannot access.}

\item{vis_graph}{A visibility graph constructed with the R package \code{pathroutr}, which is used
to reroute paths around barriers.}

\item{as_sf}{Logical. Return an \code{sf} points data frame (\code{TRUE}) or standard \code{crawl} prediction.}

\item{...}{Additional arguments passed to the \code{\link[foreach]{foreach}} function, e.g.,
for error handling in the loop.}
}
\description{
Uses a list of CRW fitted models and desired prediction times
to make location (and velocity) predictions for telemetered animals.
}
\details{
The R package \code{pathroutr} is necessary for use of the \code{barrier} rerouting.
it can be installed with the command
\code{install.packages('pathroutr', repos='https://jmlondon.r-universe.dev')}.
See 'https://github.com/jmlondon/pathroutr' for a description of use and constructing the
viability \code{vis_graph}.
}
\author{
Devin S. Johnson
}

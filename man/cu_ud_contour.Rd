% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cu_ud_contour.R
\name{cu_ud_contour}
\alias{cu_ud_contour}
\title{Approximate contours for a UD}
\usage{
cu_ud_contour(ud, prob, barrier, smoothness, ...)
}
\arguments{
\item{ud}{A \code{ud_df} object output from \code{\link{cu_kde_ud}} or  \code{\link{cu_kde_ud_sample}}.}

\item{prob}{a vector of probabilities for the contours}

\item{barrier}{A barrier polygon (e.g., land in a marine setting) that is trimmed out of the contours}

\item{smoothness}{Determines how much smoothing is done on the ud pixels.
See \code{smoothr} package function \code{\link[smoothr]{smooth}}.}

\item{...}{Further arguments (other than "smoothness") passed to \code{\link[smoothr]{smooth}}.}
}
\description{
Approximate contours for a UD
}
\author{
Devin S. Johnson
}

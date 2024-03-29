% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cu_avg_ud.R
\name{cu_avg_ud}
\alias{cu_avg_ud}
\title{Averaging Utilization Distributions}
\usage{
cu_avg_ud(ud_list, fac = NULL, w = NULL, normwt = TRUE)
}
\arguments{
\item{ud_list}{A list of individual utilization distributions calculated via
\code{\link[crawlUtils]{cu_kde_ud}} or \code{\link[crawlUtils]{cu_kde_ud_sample}}.
Each element of the list must have been calculated from the same grid created with
\code{\link[crawlUtils]{cu_ud_grid}}.}

\item{fac}{A factor variable. Averaging will be calculated for each level of \code{fac}.}

\item{w}{Weights for averaging.}

\item{normwt}{Logical. Should the weights (\code{w}) be normalized such that \code{sum(w)=1}? Defaults
to \code{TRUE}.}
}
\description{
Averaging Utilization Distributions
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crawl_utils_funcs.R
\name{cu_get_bouts}
\alias{cu_get_bouts}
\title{Label telemetry bout segments}
\usage{
cu_get_bouts(x, gap = 7, time_unit = "days")
}
\arguments{
\item{x}{column of POSIX times for which locations are observed.}

\item{gap}{specified maximum time interval for which a new bout is started. Defaults
to \code{gap=7} days.}

\item{time_unit}{Unit of time of the gap specification. Defaults to \code{"day"}.}
}
\description{
Label segments of telemetry data such that gaps between observed
locations are less than or equal to the interval specified.
}
\author{
Devin S. Johnson
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cu_crw_argos.R
\name{cu_crw_argos}
\alias{cu_crw_argos}
\title{Batch Fitting CTCRW Models for Argos (and FastGPS) Data}
\usage{
cu_crw_argos(
  data,
  move_phase = NULL,
  bm = FALSE,
  use_prior = TRUE,
  crw_control = NULL,
  fit = TRUE,
  skip_check = FALSE,
  ...
)
}
\arguments{
\item{data}{A data set containg the telemetry data for a single individual.
Will also accept a single \code{sf} data frame as well.}

\item{move_phase}{An optional character value indicating a factor variable column in the data that
designates different movement phases.}

\item{bm}{Fit a Brownian Motion model rather than in integrated OU model. Defaults to \code{bm = FALSE}.}

\item{use_prior}{Logical. Should a sensible mixture normal prior be use for the log beta and
log error scale parameters to impose a soft constraint for better numerical optimization.
Default is \code{TRUE}.}

\item{crw_control}{A named list passed to \code{\link[crawl]{crwMLE}} for optimization. Alternatives
for the default values of \code{initialSANN}, \code{attempts}, \code{control}, \code{theta}, \code{fixPar}, \code{constr}, and \code{prior}
can be specified here. See \code{\link[crawl]{crwMLE}} for a description of these arguments. WARNING!!! No
checks are made for validity of the user override. So know what you are doing.}

\item{fit}{Logical. CTCRW parameters are estimated if \code{fit=TRUE} (default), else
the results of \code{\link[crawl]{displayPar}}.}

\item{skip_check}{See \code{\link[crawl]{crwMLE}} v2.3.0. Currently ignored.}

\item{...}{Additional arguments passed to the \code{\link[foreach]{foreach}} function, e.g.,
for error handling in the loop.}
}
\description{
A basic CTCRW model is fitted to a list of data sets where each
element in the list represents the data for a single individual or deployment.
}

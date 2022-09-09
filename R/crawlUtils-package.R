#' @title Functions To Increase Usability Of The \code{crawl} Package
#'
#' @description This package is a collection of functions that enhance the \code{crawl} package for
#' for analysis of animal telemetry data. The functions integrate \code{crawl} output and the \code{sf} package
#' for ease of model fitting and track prediction, notably in marine environments.
#'
#' \tabular{ll}{
#' Package: \tab crawlUtils\cr
#' Type: \tab Package\cr
#' Version: \tab 0.1.16\cr
#' Date: \tab September 9, 2022\cr
#' License: \tab CC0 \cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' @note This software package is developed and maintained by scientists at the
#' NOAA Fisheries Pacific Islands Fisheries Science Center and should be
#' considered a fundamental research communication. The recommendations and
#' conclusions presented here are those of the authors and this software should
#' not be construed as official communication by NMFS, NOAA, or the U.S. Dept.
#' of Commerce. In addition, reference to trade names does not imply endorsement
#' by the National Marine Fisheries Service, NOAA. While the best efforts have
#' been made to insure the highest quality, tools such as this are under
#' constant development and are subject to change.
#'
#' @name crawlUtils-package
#' @aliases crawlUtils-package crawlUtils
#' @docType package
#' @author Devin S. Johnson
#' Maintainer: Devin S. Johnson <devin.johnson@@noaa.gov>
#'
NULL



.onAttach <- function(library, pkgname)
{
  info <-utils::packageDescription(pkgname)
  package <- info$Package
  version <- info$Version
  date <- info$Date
  packageStartupMessage(
    paste(package, version, paste("(",date, ")", sep=""))
  )
}

# .onUnload <- function(libpath)
# {
#   #library.dynam.unload("crawl", libpath)
#   cat("\nBye-Bye from crawl\n\n")
#   return(invisible())
# }


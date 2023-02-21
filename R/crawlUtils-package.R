#' @title Functions To Increase Usability Of The \code{crawl} Package
#'
#' @description This package is a collection of functions that enhance the \code{crawl} package for
#' for analysis of animal telemetry data. The functions integrate \code{crawl} output and the \code{sf} package
#' for ease of model fitting and track prediction, notably in marine environments.
#'
#' \tabular{ll}{
#' Package: \tab crawlUtils\cr
#' Type: \tab Package\cr
#' Version: \tab 0.1.30\cr
#' Date: \tab February 21, 2023\cr
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

rename_geometry <- function(g, name){
  current = attr(g, "sf_column")
  names(g)[names(g)==current] = name
  st_geometry(g)=name
  g
}


# weighted.var <- function(x, w, na.rm = FALSE) {
#   if (na.rm) {
#     w <- w[i <- !is.na(x)]
#     x <- x[i]
#   }
#   sum.w <- sum(w)
#   sum.w2 <- sum(w^2)
#   mean.w <- sum(x * w) / sum(w)
#   (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
#                                        na.rm)
# }

#~~~ Reference for summarizing overlapping polys:
# https://stackoverflow.com/questions/48279545/summarise-attributes-from-sfst-intersection-where-geometries-overlaps
# *** This doesn't really work in practice. Just use a common grid. ***







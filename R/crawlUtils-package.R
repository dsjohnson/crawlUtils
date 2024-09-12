#' @title Functions To Increase Usability Of The \code{crawl} Package
#'
#' @description This package is a collection of functions that enhance the \code{crawl} package for
#' for analysis of animal telemetry data. The functions integrate \code{crawl} output and the \code{sf} package
#' for ease of model fitting and track prediction, notably in marine environments.
#'
#' \tabular{ll}{
#' Package: \tab crawlUtils\cr
#' Type: \tab Package\cr
#' Version: \tab 0.1.62\cr
#' Date: \tab September 12, 2024\cr
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
#' @author Devin S. Johnson
#' Maintainer: Devin S. Johnson <devin.johnson@@noaa.gov>
# #' @useDynLib crawlUtils, .registration = TRUE
#'
"_PACKAGE"



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

#' @import dplyr
rm_dup <- function(x){
  deploy_id <- datetime <- quality <- NULL
  x <- x |>
    group_by(deploy_id) |>
    arrange(datetime, quality) |>
    mutate(
      rank = 1L,
      rank = case_when(duplicated(datetime, fromLast = FALSE) ~
                         lag(rank) + 1L, TRUE ~ rank)) |>
    dplyr::filter(rank == 1) |>
    ungroup() |>
    arrange(deploy_id, datetime)
}

#' @title Silverman default bandwidth calculation
#' @param xy Data coordinates
#' @param ess An effective sample size. If left as `NULL`, `nrow(xy)` is used.
#' @importFrom stats quantile
#' @export
bw_silver <- function (xy, ess=NULL)
{
  bw <- NULL
  if(is.null(ess)) ess <- nrow(xy)
  for(i in 1:2){
    iqr <- diff(quantile(xy[,i], c(0.25, 0.75)))
    h <- 1.06 * min(sqrt(var(xy[,i])), iqr/1.34) * ess^(-1/5)
    bw <- c(bw,h)
  }
  return(bw)
}

bw_scott <- function (xy, ess=NULL)
{
  v <- mean(var(xy[,1]), var(xy[,2]))
  if(is.null(ess)) ess <- nrow(xy)
  h <- 1.06 * sqrt(v) * ess^(-1/5)
  return(c(h,h))
}

bw_def <- function(xy, ess=NULL, groups=NULL){
  if(!is.null(groups)){
    grp_lvl <- unique(groups)
    ng <- length(grp_lvl)
    if(ng==1){
     h <- bw_scott(xy,ess)
    } else{
      h_k <- rep(NA, ng)
      for(k in 1:length(grp_lvl)){
        xy_k <- xy[groups==grp_lvl[k],]
        h_k[k] <- bw_scott(xy_k,ess)
      }
      h <- mean(h_k)
    }
  } else{
    h <- bw_scott(xy,ess)
  }
  return(h)
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







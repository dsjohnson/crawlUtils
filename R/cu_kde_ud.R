#' @title Kernel density estimator
#' @description Calculates a Gaussian KDE for \code{crawl} simulations and predictions.
#' @param pts A \code{\link[crawl]{crwPredict}} or \code{\link[crawl]{crwPostIS}} object, or
#' their 'sf' versions (See \link[crawl]{crw_as_sf}).
#' @param grid An \code{\link[sf]{sf}} data frame containing the desired grid location for UD estimation.
#' @param kern Type of covariance matrix for the Gaussian kernels.
#' @param ess Effective sample size.
#' @param norm Logical. Should each individual kernel be normalized to
#' sum-to-one over the locations in \code{grid}. Defaults to \code{kern = TRUE}
#' @param B Kernel covariance matrix. Defaults to \code{B = NULL} and a effective
#' sample size calculation is used for a plugin 2d Gaussian kernel.
#' @param type Form of the return type. Can be \code{"original"} to have the UD returned in the same form as the \code{grid}
#' argument. Or set \code{type="vector"} to return only a vector of UD values.
#' @author Devin S. Johnson
#' @import sf crawl
#' @useDynLib crawlUtils, .registration = TRUE
#' @export
#'
cu_kde_ud <- function(pts, grid, kern="iso", ess=NULL, norm=TRUE, B=NULL, type="original"){
  gorig <- NULL
  grid_id <- attr(grid, "grid_id")
  if(inherits(grid,"sf")){
    gorig <- grid
    grid <- st_geometry(gorig)
  }
  if(inherits(grid,"sfc_POLYGON")){
    grid <- sf::st_centroid(grid)
  }
  if(!inherits(grid,"sfc_POINT")){
    stop("The 'grid' argument must be either 'sfc_POLYGON', 'sfc_POINT', or 'sf' data frame containing the previous geometry types.")
  }
  if(!kern%in%c("iso","diag","full")) stop("The 'kern' argument must be one of 'iso','diag',or 'full'")
  if(inherits(pts,"crwPredict") | inherits(pts,"crwIS")){
    pts <- crw_as_sf(pts, "POINT")
  }
  if(inherits(pts,"sf")){
    if(! "TimeNum"%in%colnames(pts)) stop("It appears that 'pts' is not the correct class.")
    if(inherits(pts,"sfc_LINESTRING")) stop("The locations in 'pts' must be of class 'sfc_POINT'")
  }
  if(st_crs(pts) != st_crs(grid)) stop("The 'pts' and 'grid' crs specifications do not match.")

  if(type!="skeleton"){
    xy_grid <- st_coordinates(grid)
    xy_pts <- st_coordinates(pts)
    if(is.null(B)){
      if(kern == 'iso'){
        B <- var(c(xy_pts[,1]-mean(xy_pts[,1]), xy_pts[,2]-mean(xy_pts[,2])))*diag(2)
      } else if(kern=="diag"){
        B <- diag(diag(var(xy_pts)))
      } else{
        B <- var(xy_pts)
      }
      B <- (ess^(-1/3))*B
    }
    ud <- kde_estimate(grid=xy_grid, points=xy_pts, B=solve(B), norm = norm)
  } else {
    ud <- NA
    B <- NA
  }

  if(!is.null(ess)){
    if(inherits(ess, "crwFit")){
      ess <- cu_crw_ess(ess, pts)
    } else if(!is.numeric(ess)){
      stop("The 'ess' argument must be either a 'crwFit' object or numeric if specified.")
    }
  }else{
    ess <- nrow(pts)
  }


  if(type%in%c("original","skeleton")){
    if(!is.null(gorig)){
      gorig$ud <- ud
    } else{
      gorig <- st_as_sf(grid)
      gorig$ud <- ud
    }
    attr(gorig,"is_ud") <- TRUE
    attr(gorig, "ess") <- ess
    attr(gorig, "B") <- B
    attr(gorig, "grid_id") <- grid_id
    return(gorig)
  } else{
    attr(ud,"is_ud") <- TRUE
    attr(ud, "ess") <- ess
    attr(ud, "B") <- B
    attr(ud, "grid_id") <- grid_id
    return(as.vector(ud))
  }
}

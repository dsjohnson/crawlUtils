#' @title Kernel density estimator
#' @description Calculates a Gaussian KDE for \code{crawl} simulations and predictions.
#' @param pts A \code{\link[crawl]{crwPredict}} or \code{\link[crawl]{crwPostIS}} object, or
#' their 'sf' versions (See \link[crawl]{crw_as_sf}).
#' @param grid An \code{\link[sf]{sf}} data frame containing the desired grid location for UD estimation.
#' @param ess Effective sample size object or a `crwFit` object. See `link[crawlUtils]{cu_crw_ess}`.
#' @param use_w Use weights from the `ess` object for a weighted KDE.
#' @param norm Logical. Should each individual kernel be normalized to
#' sum-to-one over the locations in \code{grid}. Defaults to \code{kern = TRUE}
#' @param bw Kernel bandwidth (standard deviation of Gaussian kernel). Defaults to the
#' default plugin bandwidth `bandwidth.nrd` in the `MASS` package.
#' @param bw_subset A vector of values indicating which `pts` should be used for calculating `B` if left unspecified.
#' @param type Form of the return type. Can be \code{"original"} to have the UD returned in the same form as the \code{grid}
#' argument. Or set \code{type="vector"} to return only a vector of UD values.
#' @param ... additional arguments passed to \code{\link[crawlUtils]{cu_vel_B}}
#' @author Devin S. Johnson
#' @import sf crawl
#' @useDynLib crawlUtils, .registration = TRUE
#' @export
#'
cu_kde_ud <- function(pts, grid, ess=NULL, use_w=TRUE, norm=TRUE, bw=NULL, bw_subset=TRUE, type="original", ...){

  ### Checks
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
  if(inherits(pts,"crwPredict") | inherits(pts,"crwIS")){
    pts <- crw_as_sf(pts, "POINT")
  }
  if(inherits(pts,"sf")){
    if(! "TimeNum"%in%colnames(pts)) stop("It appears that 'pts' is not the correct class.")
    if(inherits(pts,"sfc_LINESTRING")) stop("The locations in 'pts' must be of class 'sfc_POINT'")
  }
  if(st_crs(pts) != st_crs(grid)) stop("The 'pts' and 'grid' crs specifications do not match.")

  ### Get ESS value
  if(!is.null(ess)){
    if(inherits(ess, "crwFit")){
      ess <- cu_crw_ess(ess, pts)
    } else if(!inherits(ess, "crwESS")){
      stop("The 'ess' argument must be either a 'crwFit' or 'crwESS' object.")
    }
  } else{
    ess <- list(Ne = nrow(pts), w = rep(1/nrow(pts), nrow(pts)))
  }

  ### Compute bandwidth matrix
  if(type!="skeleton"){
    xy_grid <- st_coordinates(grid)
    xy_pts <- st_coordinates(pts)
    if(is.null(bw)){
      defbw <- function(x,ess)
      {
        r <- quantile(x, c(0.25, 0.75))
        h <- (r[2] - r[1])/1.34
        (1.06 * min(sqrt(var(x)), h) * ess$Ne^(-1/5))^2
      }
      xy_B <- xy_pts[bw_subset,]
      B <- diag(c(defbw(xy_B[,1], ess), defbw(xy_B[,2], ess)))
    } else if(bw=="vel_B"){
      B <- cu_vel_B(pts[bw_subset,], ess, ...)
    } else if(is.numeric(bw) & length(bw==1)){
      B <- diag(rep(bw^2,2))
    }
    else{
      stop("Error in 'bw' specification!")
    }
    if(use_w){
      w <- ess$w
    } else{
      w <- rep(1/nrow(xy_pts), nrow(xy_pts))
    }
    ud <- kde_estimate(grid=xy_grid, points=xy_pts, B=solve(B), w=w, norm = norm)
  }
  else {
    ud <- NA
    B <- NA
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

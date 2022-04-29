#' @title Calculate Gaussian Kernel Utilization Distribution
#' @description Calculates a kernel density estimate using telemetry locations
#' over a grid accounting for barriers by normalizing each kernel before addition
#' to the overall sum. The final estimate will always be normalized to sum to one
#' over the estimate grid.
#' @param x an sf or sfc object. The KDE grid will be estimated from the bounding box.
#' @param barrier An sf polygon object that defines an area where use is excluded.
#' @param norm Logical. Should each individual kernel be normalized to account for barriers
#' before addition to the total KDE. Defaults to \code{TRUE} if a barrier is specified.
#' @param ... Any arguments passed to \code{sf::st_make_grid} to create the KDE prediction grid. See \link[sf:st_make_grid]{sf::st_make_grid} for
#' description of additional arguments to make the KDE grid.
#' @import sf
#' @useDynLib crawlUtils, .registration = TRUE
#' @author Devin S. Johnson
#' @export
#'
cu_kde_ud <- function(x, barrier=NULL, norm, ...){
  grid <- sf::st_make_grid(x,...)
}


#' @title Create covariance function for a fitted CRW model object
#' @description A function is created to evaluate the covariance function of the fitted
#' CRW movement model
#' @param x A crwFit object created by a call to \code{\link[crawl]{crwMLE}}
#' @details The function returns a function to evaluate the covariance of the fitted Integrated
#' Ornstein-Ulenbeck movement model. The returned function has 3 arguments: (1) \code{t1}
#' and (2) \code{t2} both vectors of times to evaluate the covariance function of the fitted
#' IOU model, and (3) \code{E} (defaults to \code{E=0}). Which is the "zero" time of the process.
#' Typically \code{E} will be the time of the first observation.
#' @references Taylor, J. M., Cumberland, W. G., & Sy, J. P. (1994).
#' A stochastic model for analysis of longitudinal AIDS data. Journal of the
#' American Statistical Association, 89(427), 727-736.
#'
#' Johnson, D. S., London, J. M., Lea, M. A., & Durban, J. W. (2008).
#' Continuous‐time correlated random walk model for animal telemetry data.
#' Ecology, 89(5), 1208-1215.
#'@export
#'@importFrom utils tail
#'@importFrom stats var
#'@author Devin S. Johnson
#'
cu_crw_covfun <- function(x){
  if((x$mov.model != ~1) | (x$random.drift==TRUE) | (!is.null(x$activity))) stop("Sorry, currently this function only works with the base CRW model")
  par <- tail(x$par,2)
  b <- exp(par[2])
  sig2 <- exp(2*par[1])
  Pvec <- diag(var(x$data[,x$coord]))
  foo <- function(t1, t2, E=min(x$data$TimeNum)){
    if(E>min(c(t1,t2))) stop("E must be < min(c(t1,t2))")
    s <- pmin(t1-E, t2-E)
    t <- pmax(t1-E, t2-E)
    sig2*(s - (1/(2*b))*(1+exp(-b*(t-s))) + (1/(2*b))*exp(-b*t) + (1/(2*b))*exp(-b*s) )
  }
  attr(foo, "b") <- b
  attr(foo, "sig2") <- sig2
  return(foo)
}

#' @title Calculate correlation matrix for a set of times from a CRW covariance function
#' @description Using a correlation function created by \code{\link{cu_crw_covfun}}
#' from a fitted CRW model a covariance (correlation) matrix is created for observations
#' at the user provided times.
#' @param x Either a \code{crwFit} object from a call to \code{\link[crawl]{crwMLE}}, or a
#' vector of times.
#' @param corr Should the function return a correlation or covariance matrix? Defaults
#' to \code{corr = TRUE}.
#' @param cf Covariance function created from a \code{crwFit} object with a call to
#' \code{\link{cu_crw_covfun}}.
#' @param E The 'zero' time used for the covariance function. Defaults to \code{E = 0}.
#' (See \code{\link{cu_crw_covfun}}).
#' @details If \code{x} is a \code{crwFit} object, then the \code{cf} and \code{E}
#' arguments are ignored. The resulting matrix is the covariance (correlation) matrix
#' for the observed location times conditioned on the first location. Therefore,
#' for n observed locations an (n-1) by (n-1) matrix will result. This is most
#' useful for the effective sample size computation for a kernel density estimate.
#' @author Devin S. Johnson
#' @importFrom stats cov2cor
#' @export
cu_crw_covmat <- function(x, corr=TRUE, cf, E=0){
  if(inherits(x, "crwFit")){
    foo <- cu_crw_covfun(x)
    E <- min(x$data$TimeNum)
    S <- outer(x$data$TimeNum[-1],x$data$TimeNum[-1], foo, E=E)
  } else {
    S <- outer(x,x,FUN=cf,E=E)
  }
  if(corr){
    S <- cov2cor(S)
  }
  return(S)
}


#' @title Calculate Effective Sample Size for a Set of CRW locations
#' @description Estimates the number of independent locations in a CRW data set
#' using the method of Acosta and Vallejos (2018) [AV18].
#' @param fit A \code{crwFit} object (See \code{\link[crawl]{crwMLE}}).
#' @param aug Either a \code{\link[crawl]{crwPredict}} or \code{\link[crawl]{crwPostIS}} objects
#' from which the extra \code{predTime} location times will be used in the calculation.
#' The \code{\link[crawl]{crw_as_sf}} transformed versions of these objects will also work.
#' @details The AV18 method was designed for spatial regression analysis, but
#' the derivations only use a general correlation matrix. Therefore, the time-series
#' correlation matrix of the CRW (IOU) process was substituted. However, there is one
#' change. The CRW (IOU) model is not stationary, so the correlation matrix of the observations
#' conditioned on the first observations is used. The resulting sample size is then
#' incremented by 1 to account for the first observation. The first observation
#' can be regarded as a single independent observation from the animal UD, the AV18
#' calculation then adds the number of additional 'independent' observations given the
#' realization of the first one.
#' @references Acosta, J., & Vallejos, R. (2018). Effective sample size for
#' spatial regression models. Electronic Journal of Statistics, 12:3147-3180.
#' @author Devin S. Johnson
#' @export
#'
cu_crw_ess <- function(fit, aug=NULL){
  df <- fit$data
  v <- diag(var(df[,fit$coord]))
  df <- fit$data[!duplicated(fit$data$TimeNum),]
  if(is.null(aug)){
    n <- nrow(df)
    R <- cu_crw_covmat(fit, corr=FALSE)
    Sx <- matrix(v[1],n,n)
    Sx[2:n,2:n] <- Sx[2:n,2:n] + R
    Sx <- cov2cor(Sx)
    Sy <- matrix(v[2],n,n)
    Sy[2:n,2:n] <- Sy[2:n,2:n] + R
    Sy <- cov2cor(Sy)
  } else{
    if(!"TimeNum"%in%colnames(aug)) stop("The 'aug' argument is not the correct class. See ?cu_crw_ess.")
    cf <- cu_crw_covfun(fit)
    first_obs <- df$TimeNum[1]==aug$TimeNum[1]
    if(!first_obs){
      times <- c(df$TimeNum[1],aug$TimeNum)
    } else{
      times <- aug$TimeNum
    }
    n <- length(times)
    R <- cu_crw_covmat(times[-1], corr=FALSE, cf=cf, E=times[1])

    Sx <- matrix(v[1],n,n)
    Sx[2:n,2:n] <- Sx[2:n,2:n] + R
    Sy <- matrix(v[2],n,n)
    Sy[2:n,2:n] <- Sy[2:n,2:n] + R
    if(!first_obs){
      Sx <- Sx[-1,-1]
      Sy <- Sy[-1,-1]
    }
    Sx <- cov2cor(Sx)
    Sy <- cov2cor(Sy)
  }
  ess <- (sum(solve(Sx,rep(1,nrow(Sx)))) + sum(solve(Sy,rep(1,nrow(Sy)))))/2
  return(ess)
}

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
#' @author Devin S. Johnson
#' @import sf crawl
#' @export
#'
cu_kern_ud <- function(pts, grid, kern="iso", ess=NULL, norm=TRUE, B=NULL){
  gorig <- NULL
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

  xy_grid <- st_coordinates(grid)
  xy_pts <- st_coordinates(pts)
  if(!is.null(ess)){
    if(inherits(ess, "crwFit")){
      ess <- cu_crw_ess(ess, pts)
    } else if(!is.numeric(ess)){
      stop("The 'ess' argument must be either a 'crwFit' object or numeric if specified.")
    }
  }else{
    ess <- nrow(pts)
  }
  if(is.null(B)){
    if(kern == 'iso'){
      B <- var(c(xy_pts[,1]-mean(xy_pts[,1]), xy_pts[,1]-mean(xy_pts[,1])))*diag(2)
    } else if(kern=="diag"){
      B <- diag(diag(var(xy_pts)))
    } else{
      B <- var(xy_pts)
    }
    B <- (ess^(-1/3))*B
  }
  ud <- kde_estimate(grid=xy_grid, points=xy_pts, B=solve(B), norm = norm)
  if(!is.null(gorig)){
    gorig$ud <- ud
  } else{
    gorig <- st_as_sf(grid)
    gorig$ud <- ud
  }
  return(gorig)
}


#' @title Create a grid for kernel density estimation of animal UDs
#' @description A \code{sf} data frame is created with polygon grid cells to evaluate a Gaussian
#' KDE at the centroids
#' @param bb An \code{sf} bounding box or something that can be coerced into a
#' bounding box via \code{\link[sf]{st_bbox}} function.
#' @param barrier \code{sf} or \code{sfc} polygon data which represents areas for
#' which the animal cannot travel. Use will only be calculated on the outside of
#' each polygon in \code{barrier}.
#' @param ... Additional arguments passed to \code{\link[sf]{st_make_grid}} which
#' is used to construct the base grid.
#' @author Devin S. Johnson
#' @import sf
#' @importFrom nngeo st_remove_holes
#' @export
#'
cu_ud_grid <- function(bb, barrier=NULL,...){
  geom <- NULL
  grid <- st_make_grid(bb, ...) %>% st_as_sf() %>%
    st_difference(barrier) %>% nngeo::st_remove_holes()
  mgrid <- filter(grid, st_is(geom,"MULTIPOLYGON")) %>% st_cast("POLYGON")
  grid <- filter(grid, st_is(geom,"POLYGON")) %>% bind_rows(mgrid)
  grid$cell <- 1:nrow(grid)
  grid$area <- st_area(grid)
  return(grid)
}










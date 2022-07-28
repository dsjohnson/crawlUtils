

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
    times <-
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
    cf <- cu_crw_covfun(fit)
    times <- df$TimeNum
    n <- length(times)
    R <- cu_crw_covmat(times[-1], corr=FALSE, cf=cf, E=times[1])
    Sx <- matrix(v[1],n,n)
    Sx[2:n,2:n] <- Sx[2:n,2:n] + R
    Sy <- matrix(v[2],n,n)
    Sy[2:n,2:n] <- Sy[2:n,2:n] + R
    Sx <- cov2cor(Sx)
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
    times <- times[!duplicated(times)]
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
#' @useDynLib crawlUtils, .registration = TRUE
#' @export
#'
cu_kde_ud <- function(pts, grid, kern="iso", ess=NULL, norm=TRUE, B=NULL){
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
  attr(gorig,"is_ud") <- TRUE
  attr(gorig, "ess") <- ess
  attr(gorig, "B") <- B
  return(gorig)
}

#' @title Kernel density estimation for a posterior sample list
#' @param smp_list A list of posterior track samples from \code{\link[crawl]{crwPostIS}} after
#' transformation via \code{\link[crawl]{crw_as_sf}}.
#' @param grid An \code{\link[sf]{sf}} data frame containing the desired grid location for UD estimation.
#' @param kern Type of covariance matrix for the Gaussian kernels.
#' @param ess Effective sample size.
#' @param norm Logical. Should each individual kernel be normalized to
#' sum-to-one over the locations in \code{grid}. Defaults to \code{kern = TRUE}
#' @param B Kernel covariance matrix. Defaults to \code{B = NULL} and a effective
#' sample size calculation is used for a plugin 2d Gaussian kernel.
#' @author Devin S. Johnson
#' @import dplyr sf
#' @importFrom stats sd
#' @export
#'
cu_kde_ud_sample <- function(smp_list, grid, kern="iso", ess=NULL, norm=TRUE, B=NULL){
  cell <- ud <- ud_tmp <- NULL
  ulist <- lapply(smp_list, cu_kde_ud, grid=grid, kern=kern, ess=ess, norm=norm, B=B)
  geom <- st_geometry(ulist[[1]])
  ulist <- lapply(ulist, st_drop_geometry)
  umat <- sapply(ulist, "[", ,"ud")
  out <- ulist[[1]]

  out$ud <- rowMeans(umat)
  out$se_ud <- apply(umat, 1, sd)
  out <- cbind(geom, out) %>% st_as_sf()
  attr(out, "is_ud") <- TRUE
  attr(out, "is_ud_smp") <- TRUE
  attr(out, "grid_id") <- attr(grid, "grid_id")
  attr(out, "ess") <- attr(ulist[[1]], "ess")
  attr(out, "B") <- attr(ulist[[1]], "B")
  return(out)
}

#' @title Averaging Utilization Distributions
#' @param ud_list A list of individual utilization distributions calculated via
#' \code{\link[crawlUtils]{cu_kde_ud}} or \code{\link[crawlUtils]{cu_kde_ud_sample}}.
#' Each element of the list must have been calculated from the same grid created with
#' \code{\link[crawlUtils]{cu_ud_grid}}.
#' @param fac A factor variable. Averaging will be calculated for each level of \code{fac}.
#' @param w Weights for averaging.
#' @import foreach sf
#' @importFrom stats weighted.mean
#' @export
cu_avg_ud <- function(ud_list, fac=NULL, w=NULL){
  i <- NULL
  if(!is.null(fac) & length(fac)!=length(ud_list)) stop("The 'fac' variable is not the same length as 'ud_list'")
  if(w=="ess"){
    w <- sapply(ud_list, attr, which="ess")
  }
  if(!is.null(w) & length(w)!=length(ud_list)) stop("The weights vector, 'w', is not the same length as 'ud_list'")
  if(is.null(w)) w <- rep(1, length(ud_list))
  if(any(w<0)) stop("There are w<0. All must be positive.")
  if(is.null(fac)) fac <- rep(1,length(ud_list))
  fac <- factor(fac)
  lfac <- levels(fac)
  ids <- sapply(ud_list, attr, which="grid_id")
  if(all(ids!=ids[1])){
    stop("UDs in 'ud_list' were not created from a 'cu_ud_grid()' grid. Averaging not possible.")
  }
  out_list <- foreach(i = 1:length(levels(fac)))%do%{
    idx <- fac==lfac[i]
    ulist <- ud_list[idx]
    geom <- st_geometry(ulist[[1]])
    ulist <- lapply(ulist, st_drop_geometry)
    umat <- sapply(ulist, "[", ,"ud")
    semat <- matrix(0, nrow(umat), ncol(umat))
    for(j in 1:ncol(umat)){
      if(!is.null(attr(ulist[[j]], "is_ud_smp"))) semat[,j] <- ulist[[j]]$se_ud
    }
    out <- ulist[[1]] %>% select(-any_of(c('ud', 'se_ud')))
    out$ud <- apply(umat, 1, weighted.mean, w=w[idx]/sum(w[idx]))
    var_ud <- apply(umat, 1, weighted.var, w=w[idx]/sum(w[idx])) + rowMeans(semat)
    out$se_ud <- sqrt(var_ud)/sqrt(length(ulist))
    out <- cbind(geom, out) %>% st_as_sf()
    attr(out, "is_ud") <- TRUE
    attr(out, "grid_id") <- ids[1]
    out
  }
  if(length(out_list)==1) out_list <- out_list[[1]]
  return(out_list)
}





weighted.var <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  sum.w2 <- sum(w^2)
  mean.w <- sum(x * w) / sum(w)
  (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
                                       na.rm)
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
#' @importFrom units set_units
#' @export
#'
cu_ud_grid <- function(bb, barrier=NULL,...){
  geom <- NULL
  if(!inherits(bb, "bbox")) bb <- st_bbox(bb)
  grid <- st_make_grid(bb, ...) %>% st_as_sf()
  if(!is.null(barrier)){
    grid <- grid %>% st_difference(barrier) %>% nngeo::st_remove_holes()
    geom <- attr(grid, "sf_column")
    mgrid <- filter(grid, st_is(.data[[geom]],"MULTIPOLYGON")) %>% st_cast("POLYGON")
    grid <- filter(grid, st_is(.data[[geom]],"POLYGON")) %>% bind_rows(mgrid)
  }
  grid$cell <- 1:nrow(grid)
  grid$area <- st_area(grid)
  if(max(grid$area)>units::set_units(1e+06,"m^2")) grid$area <- units::set_units(grid$area, "km^2")
  attr(grid,"grid_id") <- as.character(as.numeric(Sys.time()))
  return(grid)
}


#' @title Highest utilization density
#' @description The lowest 1-alpha percent of a utilization distribution is removed
#' to give a highest alpha\% UD which can be used as a home range estimate, or just
#' reduce spatial extent of UDs for a animal with a small spatial scale of use relative to the study area
#' @param ud A \code{ud_df} object output from \code{\link{cu_kde_ud}}.
#' @param alpha The percent cutoff for the highest utilization probability cells. Defaults to \code{alpha = 0.9}.
#' @author Devin S. Johnson
#' @export
cu_hud <- function(ud, alpha=0.9){
  if(is.null(attr(ud, "is_ud"))) stop("The 'ud' argument must be a UD object from the 'cu_kde_ud()' function.")
  if(zapsmall(sum(ud$ud))!=1) stop("UD values must be normalized to find HUD!")
  ud <- ud[order(ud$ud, decreasing=TRUE),]
  val <- cumsum(ud$ud)
  ud <- ud[val<=alpha,]
  ud <- ud[order(ud$cell),]
  return(ud)
}



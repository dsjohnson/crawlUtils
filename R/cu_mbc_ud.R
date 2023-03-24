#' @title Gaussian Mixture Model Utilization Density Estimation
#' @description Calculates a Gaussian KDE for \code{crawl} simulations and predictions.
#' @param pts A \code{\link[crawl]{crwPredict}} or \code{\link[crawl]{crwPostIS}} object, or
#' their 'sf' versions (See \link[crawl]{crw_as_sf}).
#' @param grid An \code{\link[sf]{sf}} data frame containing the desired grid location for UD estimation.
#' @param ess Effective sample size.
#' @param type Form of the return type. Can be \code{"original"} to have the UD returned in the same form as the \code{grid}
#' argument. Set \code{type="vector"} to return only a vector of UD values. Finally, set \code{type="gmm"}
#' to return the fitted \code{\link[mclust]{densityMclust}} object.
#' @param mclust_args A named list of additional arguments to pass to \code{\link[mclust]{densityMclust}}.
#' @author Devin S. Johnson
#' @import sf crawl mclust
#' @export
#'
cu_gmm_ud <- function(pts, grid, ess=nrow(pts), type="original", mclust_args=list()){

  ### Checks

  gorig <- NULL
  # if(debug) browser()
  if(type!="gmm"){
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
  }
  if(inherits(pts,"crwPredict") | inherits(pts,"crwIS")){
    pts <- crw_as_sf(pts, "POINT")
  }
  if(inherits(pts,"sf")){
    if(! "TimeNum"%in%colnames(pts)) stop("It appears that 'pts' is not the correct class.")
    if(inherits(pts,"sfc_LINESTRING")) stop("The locations in 'pts' must be of class 'sfc_POINT'")
  }
  if(type!="gmm"){if(st_crs(pts) != st_crs(grid)) stop("The 'pts' and 'grid' crs specifications do not match.")}

  ### Get ESS value
  if(inherits(ess, "crwFit")){
    ess <- cu_crw_ess(ess, pts)
  } else if(!is.numeric(ess)){
    stop("The 'ess' argument must be either a 'crwFit' object or numeric if specified.")
  }


  ### Compute GMM
  if(type!="skeleton"){
    ## GMM models
    # mclust_args <- as.list(match.call(expand.dots=FALSE))$...
    if(is.null(mclust_args$G)){
      G <- 1:9
    } else{
      G <- mclust_args$G
      mclust_args$G <- NULL
    }
    if(is.null(mclust_args$modelNames)){
      modelNames <- mclust.options("emModelNames")
    } else{
      modelNames <- mclust_args$modelNames
      mclust_args$modelNames <- NULL
    }
    if(is.null(mclust_args$verbose)) mclust_args$verbose <- FALSE
    if(is.null(mclust_args$plot)) mclust_args$plot <- FALSE

    likmat <- matrix(NA, nrow=length(G), ncol=length(modelNames))
    dfmat <- matrix(NA, nrow=length(G), ncol=length(modelNames))

    if(type!="gmm") xy_grid <- st_coordinates(grid)
    xy_pts <- st_coordinates(pts)
    nobs <- nrow(xy_pts)

    for(g in 1:length(G)){
      for(m in 1:length(modelNames)){
        args <-  c(list(data=xy_pts, G=G[g], modelNames=modelNames[m]), mclust_args)
        fit <- do.call(densityMclust, args)
        likmat[g,m] <- ifelse(!is.null(fit$loglik), fit$loglik, NA)*(ess/nobs)
        dfmat[g,m] <- ifelse(!is.null(fit$loglik), fit$df, NA)
      }
    }
    bic <-  -(-2*likmat + log(ess)*dfmat)
    g <- row(bic)[which.max(bic)]
    m <- col(bic)[which.max(bic)]
    args <-  c(list(data=xy_pts, G=G[g], modelNames=modelNames[m]), mclust_args)

    fit <- do.call(densityMclust, args)
    # browser()
    if(type!="gmm") {ud <- predict(fit, newdata=xy_grid, "dens"); ud <- ud/sum(ud)}
  } else {
    ud <- NA
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
    attr(gorig, "grid_id") <- grid_id
    return(gorig)
  } else if(type=="gmm"){
    # attr(fit,"is_ud") <- TRUE
    # attr(fit, "ess") <- ess
    return(fit)
  } else{
    attr(ud,"is_ud") <- TRUE
    attr(ud, "ess") <- ess
    attr(ud, "grid_id") <- grid_id
    return(as.vector(ud))
  }
}

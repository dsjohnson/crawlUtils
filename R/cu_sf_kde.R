
#' @title Spatial kernel density estimate
#' @description A weighted or unweighted Gaussian Kernel Density estimate
#'              for point spatial data
#'
#' @param x             sf POINT object
#' @param w             Optional values, associated with `x` coordinates,
#'                      to be used as weights
#' @param bw            Standard deviation scale bandwidth of Gaussian Kernel, must be units
#'                      of `x` projection.
#' @param ref           A terra SpatRaster
#' @param ess           A effective sample size to use instead of `nrow(x)` for determining the default bandwidth.
#' @param mask          (TRUE/FALSE) mask resulting raster if ref is provided
#'                      as a SpatRaster
#'
#' @details
#' Please note that ks methods for estimation has been reverted to the Gussian method proposed
#' in Venables & Ripley (2002). There was not enought evendence that the Chacon & Duong (2018)
#' multivariate method(s) for bandwidth selection and kernal estimation were suitable for
#' spatial random fields.
#'
#' @return  a terra SpatRaster class object containing kernel density estimate
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org> and Devin S. Johnson <devin.johnson@@noaa.gov>
#'
#' @references
#' Duong, T. & Hazelton, M.L. (2005) Cross-validation bandwidth matrices for multivariate
#'   kernel density estimation. Scandinavian Journal of Statistics, 32, 485-506.
#'
#' Wand, M.P. & Jones, M.C. (1994) Multivariate plug-in bandwidth selection. Computational
#'   Statistics, 9, 97-116.
#'
#' Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S.
#'   Fourth edition. Springer.
#'
#' @export
#' @import sf
#' @importFrom terra ext crds crs ncell values
#' @importFrom MASS kde2d bandwidth.nrd

cu_sf_kde <- function (x, w = NULL, bw = NULL, ref, ess=NULL, mask = FALSE)
{
  if (!inherits(x, c("sf", "sfc"))) stop(deparse(substitute(x)), " must be a sf, or sfc object")
  if (unique(as.character(sf::st_geometry_type(x))) != "POINT") stop(deparse(substitute(x)), " must be single-part POINT geometry")
  if (!inherits(ref, "SpatRaster")) stop(deparse(substitute(ref)), " must be a terra SpatRast object")

  fhat <- function(x, y, h, w, n = 25, lims = c(range(x), range(y))) {
    nx <- length(x)
    if (length(y) != nx)
      stop("data vectors must be the same length")
    if (length(w) != nx & length(w) != 1)
      stop("weight vectors must be 1 or length of data")
    if (missing(h)) {
      h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y))
    }
    else {
      h <- rep(h, length.out = 2L)
    }
    if (any(h <= 0))
      stop("bandwidths must be strictly positive")
    if (missing(w)) {
      w <- numeric(nx) + 1
    }
    gx <- seq(lims[1], lims[2], length = n[1])
    gy <- seq(lims[3], lims[4], length = n[2])
    h <- h/4
    ax <- outer(gx, x, "-")/h[1]
    ay <- outer(gy, y, "-")/h[2]
    z <- (matrix(rep(w, n[1]), nrow = n[1], ncol = nx, byrow = TRUE) *
            matrix(stats::dnorm(ax), n[1], nx)) %*% t(matrix(stats::dnorm(ay),
                                                             n[2], nx))/(sum(w) * h[1] * h[2])
    return(list(x = gx, y = gy, z = z))
  }
  if (is.null(bw)) {
    # bw <- bw_silver(st_coordinates(x), ess)
    bw <- bw_silver(st_coordinates(x), ess)
  }
  else {
    if(length(bw)==2) bw <- c(bw, bw)
  }
  ns <- c(terra::ncol(ref), terra::nrow(ref))
  xy_ref <- crds(ref, df=TRUE, na.rm=FALSE)
  x_ref <- range(xy_ref$x)
  y_ref <- range(xy_ref$y)
  if (!is.null(w)) {
    # message("\n", "calculating weighted kde", "\n")
    k <- fhat(sf::st_coordinates(x)[, 1], sf::st_coordinates(x)[,2], w = w, h = 4*bw, n = ns, lims = as.vector(c(x_ref, y_ref)))
  }
  else {
    # message("\n", "calculating unweighted kde", "\n")
    k <- MASS::kde2d(sf::st_coordinates(x)[, 1], sf::st_coordinates(x)[, 2], h = 4*bw, n = ns, lims = as.vector(c(x_ref, y_ref)))
  }
  k$z <- k$z * ncell(ref)
  zzz <- zapsmall(k$z)
  zzz <- t(zzz[, ncol(zzz):1L])
  kde.est <- terra::rast(zzz, ext=ext(ref), crs=crs(ref))
  if (mask == TRUE) {
    kde.est <- terra::mask(kde.est, ref)
  }
  # kde.est[kde.est==0] <- NA
  kde.est <- kde.est/sum(terra::values(kde.est), na.rm=TRUE)
  return(kde.est)
}

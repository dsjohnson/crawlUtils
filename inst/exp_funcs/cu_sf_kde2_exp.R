
#' @title Spatial kernel density estimate
#' @description A weighted or unweighted Gaussian Kernel Density estimate
#'              for point spatial data
#'
#' @param x             sf POINT object
#' @param w             Optional values, associated with `x` coordinates,
#'                      to be used as weights
#' @param ref           A terra SpatRaster
#' @param ess           A effective sample size to use instead of `nrow(x)` for determining the default bandwidth.
#' @param norm_kern     Logical. Should the kernels be normalized to sum to 1 over the `ref` raster before being included in the KDE. Defaults to `norm_kern=TRUE`.
#' @param mask          (TRUE/FALSE) mask resulting raster if ref is provided
#'                      as a SpatRaster
#' @param bw            Distance bandwidth of Gaussian Kernel, must be units
#'                      of projection
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
#' @importFrom terra ext crds crs ncell

cu_sf_kde <- function (x, ref, ess=NULL, norm_kern=TRUE, mask = FALSE, bw = NULL)
{
  if (!inherits(x, c("sf", "sfc"))) stop(deparse(substitute(x)), " must be a sf, or sfc object")
  if (unique(as.character(sf::st_geometry_type(x))) != "POINT") stop(deparse(substitute(x)), " must be single-part POINT geometry")
  if (!inherits(ref, "SpatRaster")) stop(deparse(substitute(ref)), " must be a terra SpatRast object")

  if (is.null(bw)) {
    h <- bw_def(st_coordinates(x), ess)
    bw <- c(h,h)
    # message("Using ", round(bw[1], 3), ", ", round(bw[2], 3), " for bandwidth", "\n")
  }
  else {
    if(length(bw)==1) bw <- c(bw, bw)
  }
  # browser()
  idx <- c(1:ncell(ref))[!is.na(terra::values(ref))]
  xy_ref <- crds(ref, na.rm=FALSE)[idx,]
  xy_pts <- st_coordinates(x)
  k <- dnorm(outer(xy_ref[,1], xy_pts[,1], "-")/bw[1]) * dnorm(outer(xy_ref[,2], xy_pts[,2], "-")/bw[2])
  if(norm_kern) k <- sweep(k, 2, colSums(k), "/")
  k <- zapsmall(rowSums(k)*ncell(ref))
  kde.est <- rast(ref)
  kde.est[idx] <- k
  kde.est <- kde.est/sum(values(kde.est), na.rm=TRUE)
  return(kde.est)
}

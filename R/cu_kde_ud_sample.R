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


#' @title Kernel density estimation for a posterior sample list
#' @param smp_list A list of posterior track samples from \code{\link[crawl]{crwPostIS}} after
#' transformation via \code{\link[crawl]{crw_as_sf}}.
#' @param grid An \code{\link[sf]{sf}} data frame containing the desired grid location for UD estimation.
#' @param ess Effective sample size.
#' @param return_gmm Logical. Return a list of Gaussian mixture model fits. Default is \code{FALSE}.
#' @param mclust_args Named list of additional arguments to pass to \code{\link[mclust]{densityMclust}}
#' @author Devin S. Johnson
#' @import dplyr sf
#' @importFrom stats sd
#' @export
#'
cu_gmm_ud_sample <- function(smp_list, grid, ess=NULL, return_gmm=FALSE, mclust_args=list()){
  if(!return_gmm){
    out <- cu_gmm_ud(smp_list[[1]], grid=grid, ess=ess, type="skeleton", mclust_args=mclust_args)
    umat <- sapply(smp_list, cu_gmm_ud, grid=grid, ess=ess, type="vector", mclust_args=mclust_args)
    out$ud <- rowMeans(umat)
    out$se_ud <- apply(umat, 1, sd)
  } else{
    out <- lapply(smp_list, cu_gmm_ud, grid=grid, ess=ess, type="gmm", mclust_args=mclust_args)
  }
  return(out)
}

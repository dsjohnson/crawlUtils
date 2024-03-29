#' @title Kernel density estimation for a posterior sample list
#' @param smp_list A list of posterior track samples from \code{\link[crawl]{crwPostIS}} after
#' transformation via \code{\link[crawl]{crw_as_sf}}.
#' @param grid An \code{\link[sf]{sf}} data frame containing the desired grid location for UD estimation.
#' @param ess Effective sample size.
#' @param norm Logical. Should each individual kernel be normalized to
#' sum-to-one over the locations in \code{grid}. Defaults to \code{kern = TRUE}
#' @param bw Kernel covariance matrix. Defaults to \code{B = NULL} and a effective
#' sample size calculation is used for a plugin 2d Gaussian kernel.
#' @param bw_subset A vector of values indicating which `pts` should be used for calculating `bw` if left unspecified.
#' @param ... Additional arguments passed to \code{\link[crawlUtils]{cu_vel_B}}
#' @author Devin S. Johnson
#' @import dplyr sf
#' @importFrom stats sd
#' @export
cu_kde_ud_sample <- function(smp_list, grid, ess=NULL, use_w=TRUE, norm=TRUE, bw=NULL, bw_subset=TRUE,...){
  cell <- ud <- ud_tmp <- NULL
  out <- cu_kde_ud(smp_list[[1]], grid=grid, ess=ess, norm=norm, bw=bw, use_w=use_w, bw_subset=bw_subset, type="skeleton",...)
  umat <- sapply(smp_list, cu_kde_ud, grid=grid, ess=ess, norm=norm, bw=bw, use_w=use_w, bw_subset=bw_subset, type="vector",...)
  out$ud <- rowMeans(umat)
  out$se_ud <- apply(umat, 1, sd)
  return(out)
}

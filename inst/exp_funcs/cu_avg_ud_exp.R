#' @title Averaging Utilization Distributions
#' @param ud_list A list of individual utilization distributions calculated via
#' \code{\link[crawlUtils]{cu_kde_ud}} or \code{\link[crawlUtils]{cu_kde_ud_sample}}.
#' Each element of the list must have been calculated from the same grid created with
#' \code{\link[crawlUtils]{cu_ud_grid}}.
#' @param fac A factor variable. Averaging will be calculated for each level of \code{fac}.
#' @param w Weights for averaging.
#' @param normwt Logical. Should the weights (\code{w}) be normalized such that \code{sum(w)=1}? Defaults
#' to \code{TRUE}.
#' @import foreach sf
#' @importFrom stats weighted.mean
#' @export
cu_avg_ud <- function(ud_list, fac=NULL, w=NULL, normwt=TRUE){
  i <- NULL
  if(!is.null(fac) & length(fac)!=length(ud_list)) stop("The 'fac' variable is not the same length as 'ud_list'")
  if(w=="ess"){
    w <- sapply(ud_list, attr, which="ess")
  }
  if(!is.null(w) & length(w)!=length(ud_list)) stop("The weights vector, 'w', is not the same length as 'ud_list'")
  if(is.null(w)) w <- rep(1, length(ud_list))
  if(any(w<0)) stop("There are w<0. All must be positive.")
  if(normwt) w <- w/sum(w)
  if(is.null(fac)) fac <- rep(1,length(ud_list))
  fac <- factor(fac)
  lfac <- levels(fac)
  ids <- sapply(ud_list, attr, which="grid_id")
  if(all(ids!=ids[1])){
    stop("UDs in 'ud_list' were not created from the same 'cu_ud_grid()' grid. See ?cu_average_ud_2()")
  }
  out_list <- foreach(i = 1:length(levels(fac)))%do%{
    idx <- fac==lfac[i]
    ulist <- ud_list[idx]
    geom <- st_geometry(ulist[[1]])
    ulist <- lapply(ulist, st_drop_geometry)
    umat <- sapply(ulist, "[", ,"ud")
    varmat <- matrix(0, nrow(umat), ncol(umat))
    for(j in 1:ncol(umat)){
      if(!is.null(attr(ulist[[j]], "is_ud_smp"))) varmat[,j] <- ulist[[j]]$se_ud^2
    }
    out <- ulist[[1]] %>% select(-any_of(c('ud', 'se_ud')))
    out$ud <- apply(umat, 1, crossprod, y=w)
    var_ud <- apply(varmat, 1, crossprod, y=w^2)
    out$se_ud <- sqrt(var_ud)
    out <- cbind(geom, out) %>% st_as_sf()
    attr(out, "is_ud") <- TRUE
    attr(out, "grid_id") <- ids[1]
    out
  }
  if(length(out_list)==1) out_list <- out_list[[1]]
  return(out_list)
}

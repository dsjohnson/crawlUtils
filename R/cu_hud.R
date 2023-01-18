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

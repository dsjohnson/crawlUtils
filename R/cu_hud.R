# #' @title Highest utilization density
# #' @description The lowest 1-prob percent of a utilization distribution is removed
# #' to give a highest prob\% UD which can be used as a home range estimate, or just
# #' reduce spatial extent of UDs for a animal with a small spatial scale of use relative to the study area
# #' @param ud A \code{ud_df} object output from \code{\link{cu_kde_ud}}.
# #' @param prob The percent cutoff for the highest utilization probability cells. Defaults to \code{prob = 0.9}.
# #' @author Devin S. Johnson
# #' @export
# cu_hud <- function(ud, prob=0.9){
#   if(is.null(attr(ud, "is_ud"))) stop("The 'ud' argument must be a UD object from the 'cu_kde_ud()' function.")
#   if(zapsmall(sum(ud$ud))!=1) stop("UD values must be normalized to find HUD!")
#   ud <- ud[order(ud$ud, decreasing=TRUE),]
#   val <- cumsum(ud$ud)
#   ud <- ud[val<=prob,]
#   ud <- ud[order(ud$cell),]
#   return(ud)
# }

#' @title Highest utilization density
#' @description The lowest 1-alpha percent of a utilization distribution is removed
#' to give a highest prob\% UD which can be used as a home range estimate, or just
#' reduce spatial extent of UDs for a animal with a small spatial scale of use relative to the study area
#' @param ud A `SpatRaster` object output from \code{\link{cu_sf_kde}}.
#' @param prob The percent cutoff for the highest utilization probability cells. Defaults to \code{prob = 0.9}.
#' @author Devin S. Johnson
#' @export
cu_hud <- function(ud, prob=0.9){
  if(!inherits(ud, 'SpatRaster')) stop("The 'ud' argument must be a 'SpatRaster' from the 'cu_sf_kde()' function.")
  ud_v <- cbind(1:ncell(ud), values(ud))
  ud_v <- ud_v[!is.na(ud_v[,2]),]
  ud_v[,2] <- zapsmall(ud_v[,2]/sum(ud_v[,2]))
  ud_v <- ud_v[order(ud_v[,2], decreasing=TRUE),]
  val <- cumsum(ud_v[,2])
  idx_na <- ud_v[val>prob,][,1]
  ud[idx_na] <- NA
  return(ud)
}

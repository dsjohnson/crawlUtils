#' @title Approximate contours for a UD
#' @param ud A \code{ud_df} object output from \code{\link{cu_kde_ud}} or  \code{\link{cu_kde_ud_sample}}.
#' @param prob a vector of probabilities for the contours
#' @param barrier A barrier polygon (e.g., land in a marine setting) that is trimmed out of the contours
#' @param smoothness Determines how much smoothing is done on the ud pixels.
#' See {smoothr} package function \code{\link[smoothr]{smooth}}.
#' @param ... Further arguments (other than "smoothness") passed to \code{\link[smoothr]{smooth}}.
#' @author Devin S. Johnson
#' @import smoothr
#' @export
cu_ud_contour <- function(ud, prob, barrier, smoothness, ...){
  if(!attr(ud, "is_ud")) stop("'ud' is not a {crawlUtils} UD object!")
  if(missing(prob)) prob <- c(0.25, 0.5, 0.75, 0.95)
  if(missing(smoothness)) smoothness <- 10
  c_df <- data.frame(prob = prob)
  geom <- vector("list", length(prob))
  for(p in 1:length(prob)){
    hud <- cu_hud(ud, prob[p]) |> st_geometry() |> st_union()
    geom[[p]] <- smooth(hud, method = "ksmooth", smoothness=smoothness,...) |> st_make_valid()
    # if(!missing(barrier)) geom[[p]] <- st_difference(geom[[p]], barrier)
    geom[[p]] <- st_cast(geom[[p]], "MULTILINESTRING")
  }
  c_df$geometry <- do.call("c", geom)
  c_df <- st_as_sf(c_df)
  return(c_df)
}

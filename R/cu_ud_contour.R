# #' @title Approximate contours for a UD
# #' @param ud A \code{ud_df} object output from \code{\link{cu_kde_ud}} or  \code{\link{cu_kde_ud_sample}}.
# #' @param prob a vector of probabilities for the contours
# #' @param smoothness Determines how much smoothing is done on the ud pixels.
# #' See \code{smoothr} package function \code{\link[smoothr]{smooth}}.
# #' @param ... Further arguments (other than "smoothness") passed to \code{\link[smoothr]{smooth}}.
# #' @author Devin S. Johnson
# #' @import smoothr
# #' @export
# cu_ud_contour <- function(ud, prob, smoothness, ...){
#   if(!attr(ud, "is_ud")) stop("'ud' is not a {crawlUtils} UD object!")
#   if(missing(prob)) prob <- c(0.25, 0.5, 0.75, 0.95)
#   if(missing(smoothness)) smoothness <- 10
#   c_df <- data.frame(prob = prob)
#   geom <- vector("list", length(prob))
#   for(p in 1:length(prob)){
#     hud <- cu_hud(ud, prob[p]) |> st_geometry() |> st_union()
#     geom[[p]] <- smooth(hud, method = "ksmooth", smoothness=smoothness,...) |> st_make_valid()
#     # if(!missing(barrier)) geom[[p]] <- st_difference(geom[[p]], barrier)
#     # geom[[p]] <- st_cast(geom[[p]], "MULTILINESTRING")
#   }
#   c_df$geometry <- do.call("c", geom)
#   c_df <- st_as_sf(c_df)
#   return(c_df)
# }

#' @title Approximate contours for a UD
#' @param ud A `SpatRaster` created by `\link[crawlUtils]{cu_sf_kde}`
#' @param prob a vector of probabilities for the contours
#' @param trim_pixel Logical. Should contour areas less than a pixel size be trimmed? defaults to `TRUE`.
#' @author Devin S. Johnson
#' @import sf
#' @importFrom terra as.contour res
#' @export
cu_ud_contour <- function(ud, prob=seq(0.9,0.1,-0.1), trim_pixel=TRUE){
  if(!inherits(ud, "SpatRaster")) stop("'ud' is not a 'SpatRaster' object!")
  c_df <- data.frame(prob = prob)
  geom <- vector("list", length(prob))
  ud[is.na(ud)] <- 0
  ud_v <- cbind(1:ncell(ud), values(ud))
  ud_v <- ud_v[!is.na(ud_v[,2]),]
  ud_v[,2] <- zapsmall(ud_v[,2]/sum(ud_v[,2]))
  ud_v <- ud_v[order(ud_v[,2], decreasing=TRUE),]
  val <- cumsum(ud_v[,2])
  for(p in 1:length(prob)){
    lvl <- min(ud_v[val<=prob[p],2])
    suppressWarnings({
      geom[[p]] <- st_as_sf(as.contour(ud, levels=lvl)) |> st_cast('LINESTRING') |>
        st_cast("POLYGON") |> st_make_valid() %>% st_geometry()
    })
    if(trim_pixel){
      a <- st_area(geom[[p]])
      geom[[p]] <- geom[[p]][a > set_units(prod(res(ud)), "m^2")]
    }
    geom[[p]] <- st_union(geom[[p]])
  }
  c_df$geometry <- do.call("c", geom)
  c_df <- st_as_sf(c_df)
  return(c_df)
}

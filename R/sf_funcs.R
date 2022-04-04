
#' @title Expand Spatial Bounding Box
#' @description Expand an \code{sf} bounding box by an expansion factor
#' @param bbox An \code{sf} bounding box. See \code{\link[sf:st_bbox]{sf::st_bbox}}.
#' @param ef Expansion factor, must be positive
#' @author Josh M. London
#' @importsFrom sf st_bbox
#' @export
#'
st_expand <- function(bbox, ef) {
  xmin <- as.numeric(bbox$xmin)
  xmax <- as.numeric(bbox$xmax)
  ymin <- as.numeric(bbox$ymin)
  ymax <- as.numeric(bbox$ymax)
  x_min <- xmin - ef*(xmax-xmin)
  x_max <- xmax + ef*(xmax-xmin)
  y_min <- ymin - ef*(ymax-ymin)
  y_max <- ymax + ef*(ymax-ymin)
  bbox <- st_bbox(c(xmin = x_min, xmax = x_max,
                    ymax = y_max, ymin = y_min),
                  crs = st_crs(bbox))
  return(bbox)
}


#' @title Predicate function for st_filter
#' @description Predicate function to use with \code{st_filter} such that
#' such that elemets of one spatial object are selected if
#' they are not contained at all in the other. See \code{\link[sf:st_within]{sf::st_within}}
#' @param x object of class sf, sfc or sfg
#' @param y object of class sf, sfc or sfg; if missing, x is used
#' @param sparse ogical; should a sparse index list be returned (TRUE) or a dense logical matrix? See \link[sf:st_within]{sf::st_within}.
#' @param prepared ogical; prepare geometry for x, before looping over y? See \link[sf:st_within]{sf::st_within}.
#' @param ... passed on to s2_options
#' @import sf
#' @export
#'
st_not_within <- function(x,y,sparse=TRUE,prepared=TRUE,...){
  !sf::st_within(x,y,sparse,prepared,...)
}

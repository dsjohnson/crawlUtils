
#' @title Expand Spatial Bounding Box
#' @description Expand an \code{sf} bounding box by an expansion factor
#' @param bbox An \code{sf} bounding box. See \code{\link[sf:st_bbox]{sf::st_bbox}}.
#' @param ef Expansion factor, must be positive
#' @author Josh M. London
#' @importFrom sf st_bbox
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

#' @title Calculate cellsize value for hexagon grid
#' @description Calculates the appropriate \code{cellsize} argument for making
#' a hexagon grid with \code{\link[sf]{st_make_grid}}.
#' @param area A value (m^2) for the resulting area of a full hexagon cell
#' @param radius The value for the distance (m) from the centroids to the edge of full hexagon cells.
#' @param sep The distance (m) between centoids of the hexagon grid.
#' @author Devin S. Johnson
#' @references See \url{https://github.com/r-spatial/sf/issues/1505}
#' @importFrom units set_units
#' @export
#'
hex_size <- function(area=NULL, radius=NULL, sep=NULL){

  if(!is.null(area)){
    area <- units::set_units(area, "m^2")
    return(as.numeric(2 * sqrt(area/((3*sqrt(3)/2))) * sqrt(3)/2))
  }
  if(!is.null(radius)){
    radius <- units::set_units(radius, "m")
    return(as.numeric(2*radius/sqrt(3)))
  }
  if(!is.null(sep)){
    sep <- units::set_units(sep, "m")
    return(as.numeric(sep/sqrt(3)))
  }
  stop("Argument not specified.")
}

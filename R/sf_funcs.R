
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


#' @title \code{sf::st_bbox} for a list of \code{sf} or \code{sfc} objects.
#'@param x A list of \code{sf} or \code{sfc} objects.
#'@param union Logical. Should the bounding box of the union be returned instead of
#'a list of bounding boxes.
#'@param as_sfc Logical. Should the bounding box (boxes) be returned as \code{sfc} objects.
#'@importFrom sf st_bbox st_as_sfc
#'@export
#'@author Devin S. Johnson
#'
st_bbox_list <- function(x, union=TRUE, as_sfc=FALSE){
  out <- lapply(x, st_bbox)
  if(union){
    out <- lapply(out, st_as_sfc)
    out <- st_union_list(out)
    out <- st_bbox(out)
    if(as_sfc) out <- st_as_sfc(out)
  } else{
    if(as_sfc) out <- lapply(out, st_as_sfc)
  }
  return(out)
}

#'@title  \code{sf::st_union} for a list of \code{sf} or \code{sfc} objects.
#'@param x A list of \code{sf} or \code{sfc} objects.
#'@importFrom sf st_union
#'@export
#'@author Devin S. Johnson
#'
st_union_list <- function(x){
  return(st_union(do.call(c, x)))
}

#' @title Predicate function for st_filter
#' @description Predicate function to use with \code{st_filter} such that
#' such that elements of one spatial object are selected if
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


#' @title Convert longitude/latitude coordinates from -180/180 to 0/360
#' @description Converts sf data with EPSG = 4326 from -180/180 specification to
#' 0/360 for plotting with the mapview package etc.
#' @param x An sf data frame with EPSG=4326.
#' @export
#' @import sf
#' @author Josh London
st_to_360 <- function(x){
  if(!st_crs(x)$epsg==4326) stop("This funtion is only applicable for data with EPSG=4326!")
  coords <- (sf::st_geometry(x) + c(360,90)) %% c(360) - c(0,90)
  x <- sf::st_set_geometry(x,coords) %>% sf::st_set_crs(4326)
  return(x)
}

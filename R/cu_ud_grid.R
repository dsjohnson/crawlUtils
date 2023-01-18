#' @title Create a grid for kernel density estimation of animal UDs
#' @description A \code{sf} data frame is created with polygon grid cells to evaluate a Gaussian
#' KDE at the centroids
#' @param bb An \code{sf} bounding box or something that can be coerced into a
#' bounding box via \code{\link[sf]{st_bbox}} function.
#' @param barrier \code{sf} or \code{sfc} polygon data which represents areas for
#' which the animal cannot travel. Use will only be calculated on the outside of
#' each polygon in \code{barrier}.
#' @param ... Additional arguments passed to \code{\link[sf]{st_make_grid}} which
#' is used to construct the base grid.
#' @author Devin S. Johnson
#' @import sf
#' @importFrom nngeo st_remove_holes
#' @importFrom units set_units
#' @export
#'
cu_ud_grid <- function(bb, barrier=NULL,...){
  geom <- NULL
  if(!inherits(bb, "bbox")) bb <- st_bbox(bb)
  grid <- st_make_grid(bb, ...) %>% st_as_sf()
  if(!is.null(barrier)){
    grid <- grid %>% st_difference(barrier) %>% nngeo::st_remove_holes()
    geom <- attr(grid, "sf_column")
    mgrid <- filter(grid, st_is(.data[[geom]],"MULTIPOLYGON")) %>% st_cast("POLYGON")
    grid <- filter(grid, st_is(.data[[geom]],"POLYGON")) %>% bind_rows(mgrid)
  }
  grid$cell <- 1:nrow(grid)
  grid$area <- st_area(grid)
  if(max(grid$area)>units::set_units(1e+06,"m^2")) grid$area <- units::set_units(grid$area, "km^2")
  attr(grid,"grid_id") <- as.character(as.numeric(Sys.time()))
  return(grid)
}

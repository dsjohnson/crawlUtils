#' @title Create a grid for kernel density estimation of animal UDs
#' @description A \code{sf} data frame is created with polygon grid cells to evaluate a Gaussian
#' KDE at the centroids
#' @param bb An \code{sf} bounding box or something that can be coerced into a
#' bounding box via \code{\link[sf]{st_bbox}} function.
#' @param barrier \code{sf} or \code{sfc} polygon data which represents areas for
#' which the animal cannot travel. Use will only be calculated on the outside of
#' each polygon in \code{barrier}.
#' @param remove_holes Remove holes in grid. Can sometimes speed up computations if small islands or lakes are removed.
#' defaults to `TRUE`.
#' @param ... Additional arguments passed to \code{\link[sf]{st_make_grid}} which
#' is used to construct the base grid.
#' @author Devin S. Johnson
#' @import sf
#' @importFrom nngeo st_remove_holes
#' @importFrom units set_units
#' @export
#'
cu_ud_grid <- function(bb, barrier=NULL, remove_holes=TRUE, ...){
  geom <- NULL
  if(!inherits(bb, "bbox")) bb <- st_bbox(bb)
  grid <- st_make_grid(bb, ...) %>% st_as_sf()

  # browser()

  if(!is.null(barrier)){
    barrier <- st_union(barrier)
    grid_bb <- st_bbox(grid) %>% st_expand(1.1) %>% st_as_sfc()
    barrier <- st_intersection(barrier, grid_bb)
    idx <- lengths(st_intersects(grid, barrier)) > 0
    grid_nc <- grid[!idx,] %>% rename_geometry("geometry")
    grid_c <- grid[idx,] %>%  rename_geometry("geometry")
    rm(grid)
    grid_c <- grid_c %>% st_difference(barrier)
    if(remove_holes) grid_c <- nngeo::st_remove_holes(grid_c)
    geom <- attr(grid_c, "sf_column")
    mgrid_c <- filter(grid_c, st_is(.data[[geom]],"MULTIPOLYGON")) %>% st_cast("POLYGON")
    grid_c <- filter(grid_c, st_is(.data[[geom]],"POLYGON")) %>% bind_rows(mgrid_c)
    rm(mgrid_c)
    grid_c <- rename_geometry(grid_c, "geometry")
    grid <- bind_rows(grid_nc, grid_c) %>% st_geometry() %>% st_as_sf()
    grid <- rename_geometry(grid, "geometry")
  }
  grid$cell <- 1:nrow(grid)
  grid$area <- st_area(grid)
  if(max(grid$area)>units::set_units(1e+06,"m^2")) grid$area <- units::set_units(grid$area, "km^2")
  attr(grid,"grid_id") <- as.character(as.numeric(Sys.time()))
  return(grid)
}

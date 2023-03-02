#' @title Create a grid for kernel density estimation of animal UDs
#' @description A \code{sf} data frame is created with polygon grid cells to evaluate a Gaussian
#' KDE at the centroids
#' @param bb An \code{sf} bounding box or something that can be coerced into a
#' bounding box via \code{\link[sf]{st_bbox}} function.
#' @param barrier \code{sf} or \code{sfc} polygon data which represents areas for
#' which the animal cannot travel. Use will only be calculated on the outside of
#' each polygon in \code{barrier}.
#' @param clique Logical. Should separated segments be noted so they might be removed later. Only applicable
#' when a barrier is provided.
#' @param ... Additional arguments passed to \code{\link[sf]{st_make_grid}} which
#' is used to construct the base grid.
#' @author Devin S. Johnson
#' @import sf
#' @importFrom nngeo st_remove_holes
#' @importFrom units set_units
#' @export
#'
cu_ud_grid <- function(bb, barrier=NULL, clique=FALSE, ...){
  geom <- NULL
  if(!inherits(bb, "bbox")) bb <- st_bbox(bb)
  grid <- st_make_grid(bb, ...) %>% st_as_sf()
  if(!is.null(barrier)){
    grid_bb <- st_bbox(grid) %>% st_expand(1.1) %>% st_as_sfc()
    barrier <- st_intersection(barrier, grid_bb)
    idx <- lengths(st_intersects(grid, barrier)) > 0
    grid_nc <- grid[!idx,] %>% rename_geometry("geometry")
    grid_c <- grid[idx,]
    rm(grid)
    grid_c <- grid_c %>% st_difference(barrier) %>% nngeo::st_remove_holes()
    geom <- attr(grid_c, "sf_column")
    mgrid_c <- filter(grid_c, st_is(.data[[geom]],"MULTIPOLYGON")) %>% st_cast("POLYGON")
    grid_c <- filter(grid_c, st_is(.data[[geom]],"POLYGON")) %>% bind_rows(mgrid_c)
    rm(mgrid_c)
    grid_c <- rename_geometry(grid_c, "geometry")
    grid <- bind_rows(grid_nc, grid_c)
  }
  grid$cell <- 1:nrow(grid)
  grid$area <- st_area(grid)
  if(!is.null(barrier) & clique){
    message("The 'clique' argument is not currnetly functional. It is ignored for the time being.")
  #   dmat <- st_distance(grid) |> units::drop_units()
  #   hc <- hclust(as.dist(dmat>1), method="single")
  #   grid$clique = cutree(hc, h=0.5)
  }
  if(max(grid$area)>units::set_units(1e+06,"m^2")) grid$area <- units::set_units(grid$area, "km^2")
  attr(grid,"grid_id") <- as.character(as.numeric(Sys.time()))
  return(grid)
}

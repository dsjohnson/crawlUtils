#' @title Get A Coastline \code{sf} Polygon Object For Plotting and Mapping
#' @description Uses downloaded OSM data for constructing an \code{sf} polygon coastline data object.
#' Prior to using this function you must run \code{crawlUtils::cu_download_osm()}.
#' See \code{\link[crawlUtils:cu_download_osm]{crawlUtils::cu_download_osm}}.
#' @param x An \code{sf} spatial object. The coastline will be cropped to the bounding box.
#' @param keep The amount of data retained after simplification with \code{\link[rmapshaper:ms_simplify]{rmapshaper::ms_simplify}}
#' @param union Logical. Should the retuned object be returned as a single \code{sf} \code{MULTIPOLYGON} object
#' @author Josh M. London and Devin S. Johnson
#' @import sf
#' @importFrom rmapshaper ms_simplify
#' @export
#'
cu_osm_coast <- function(x, keep=0.2, union=TRUE) {
  message("This can take a minute... maybe grab some coffee?")
  rds_file <- file.path(system.file(package="crawlUtils"), "inst", "osm", "land-polygons-complete-4326", "land.rds")
  shp_file <- file.path(system.file(package="crawlUtils"), "inst", "osm", "land-polygons-complete-4326", "land_polygons.shp")

  message("Reading in world polygon data ...")
  if(file.exists(rds_file)){
    land <- readRDS(rds_file)
  } else if(file.exists(shp_file)){
    land <- st_read(shp_file)
  }  else{
    stop("OSM data has not been previously downloaded.\nSee ?crawlUtils::cu_download_osm")
  }

  crs_land <- sf::st_crs(land)
  x_crs <- sf::st_crs(x)
  x <- x |> st_geometry() |> st_bbox() |> st_as_sfc() |> st_transform(st_crs(land))

  # if(sign(x_bb$xmin) != sign(x_bb$xmax)){
  #   # the following section includes steps to crop across the 180.
  #   target_crs <- st_crs("+proj=eqc +x_0=0 +y_0=0 +lat_0=0 +lon_0=175")
  #   # define a long & slim polygon that overlaps the meridian line & set its CRS to match
  #   # that of world
  #   # Centered in lon 175
  #   offset <- 180 - 175
  #   polygon <- st_polygon(x = list(rbind(
  #     c(-0.0001 - offset, 90),
  #     c(0 - offset, 90),
  #     c(0 - offset, -90),
  #     c(-0.0001 - offset, -90),
  #     c(-0.0001 - offset, 90)
  #     ))) %>% st_sfc() %>% st_set_crs(4326)
  #   # modify world dataset to remove overlapping portions with world's polygons
  #   osm_land <- osm_land %>% st_difference(polygon)
  #   # Transform
  #   osm_land <- osm_land %>% st_transform(crs = target_crs)
  #   x_bb <- x %>% st_convex_hull() %>% st_transform(target_crs) %>% st_bbox()
  #   message("180 adjusted...")
  # }

  message("Cropping world to 'x' bounding box ...")
  land <- suppressWarnings(st_intersection(land, x))
  land <- st_transform(land, crs=x_crs)

  message("Simplifying ...")
  if(keep>0 & keep<1) land <- rmapshaper::ms_simplify(land, keep=keep)

  message("Unioning poolygons ...")
  if(union) land <- st_geometry(land) %>% st_union()

  message("Done ...")
  return(land)
}



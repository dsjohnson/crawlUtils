#' @title Download Open Street Map Coastline Shapefile
#' @description Downloads OSM data for plotting telemetry data and use distributions.
#' @param force Logical. If data has previously been downloaded, it will force a new
#' download and update of the OSM shapefile.
#' @references https://osmdata.openstreetmap.de/data/land-polygons.html
#' @author Josh M. London and Devin S. Johnson
#' @importFrom utils download.file unzip
#' @export
#'
cu_download_osm <- function(force = FALSE) {
  dir <- file.path(system.file(package="crawlUtils"), "inst", "osm")
  if(file.exists(file.path(dir,"land-polygons-complete-4326", "land_polygons.shp")) & !force) {
    stop("OSM data has already been downloaded. Use 'force=TRUE' to update.")
  } else {
    inp <- readline(prompt = "This function will download a considerable amount of coastline data.\nAre you sure you want to proceed? [y/n]: ")
    if(tolower(inp)%in%c("n","no")) stop("OSM download crisis averted, phewww!")
    options(timeout = max(1000, getOption("timeout")))
    dir.create(dir, recursive=TRUE)
    download.file("https://osmdata.openstreetmap.de/download/land-polygons-complete-4326.zip",
                  destfile = file.path(dir,"land-polygons-complete-4326.zip"),
                  method = "auto")
    unzip(file.path(dir,"land-polygons-complete-4326.zip"), exdir = dir)
    file.remove(file.path(dir,"land-polygons-complete-4326.zip"))
    m <- paste(readLines(file.path(dir, "land-polygons-complete-4326","README.txt")),"\n")
    message(m)
  }
}

#' @title Get A Coastline \code{sf} Polygon Object For Plotting and Mapping
#' @description Uses downloaded OSM data for constructing an \code{sf} polygon coastline data object.
#' Prior to using this function you must run \code{crawlUtils::cu_download_osm()}.
#' See \code{\link[crawlUtils:cu_download_osm]{crawlUtils::cu_download_osm}}.
#' @param x An \code{sf} spatial object. The coastline will be cropped to the bounding box.
#' @param keep The amount of data retained after simplification with \code{\link[rmapshaper:ms_simplify]{rmapshaper::ms_simplify}}
#' @author Josh M. London and Devin S. Johnson
#' @import sf
#' @importFrom rmapshaper ms_simplify
#' @export
#'
cu_osm_coast <- function(x, keep=0.2) {
  osm_SHP_file <- file.path(system.file(package="crawlUtils"), "inst", "osm","land-polygons-complete-4326","land_polygons.shp")
  if(!file.exists(osm_SHP_file)) stop("OSM data has not been previously downloaded.\nSee ?crawlUtils::download_osm")
  message("This can take a while... maybe grab some coffee?")
  osm_land <- sf::read_sf(osm_SHP_file)
  x_bb <- x %>% st_convex_hull() %>% st_transform(4326) %>% st_bbox()
  if(sign(x_bb$xmin) != sign(x_bb$xmax)){
    # the following section includes steps to crop across the 180.
    target_crs <- st_crs("+proj=eqc +x_0=0 +y_0=0 +lat_0=0 +lon_0=175")
    # define a long & slim polygon that overlaps the meridian line & set its CRS to match
    # that of world
    # Centered in lon 175
    offset <- 180 - 175
    polygon <- st_polygon(x = list(rbind(
      c(-0.0001 - offset, 90),
      c(0 - offset, 90),
      c(0 - offset, -90),
      c(-0.0001 - offset, -90),
      c(-0.0001 - offset, 90)
      ))) %>% st_sfc() %>% st_set_crs(4326)
    # modify world dataset to remove overlapping portions with world's polygons
    osm_land <- osm_land %>% st_difference(polygon)
    # Transform
    osm_land <- osm_land %>% st_transform(crs = target_crs)
    x_bb <- x %>% st_convex_hull() %>% st_transform(target_crs) %>% st_bbox()
    message("180 adjusted...")
  }
  osm_land <- st_crop(osm_land, x_bb)
  message("cropped...")
  osm_land <- st_transform(osm_land, crs=st_crs(x))
  message("projected...")
  osm_land <- rmapshaper::ms_simplify(osm_land, keep=keep)
  message("simplified...")
  osm_land <- st_union(osm_land)
  message("unioned...Done...")
  return(osm_land)
}



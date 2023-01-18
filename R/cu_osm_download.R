#' @title Download Open Street Map Coastline Shapefile from Daylight Map Distribution
#' @description Downloads OSM data for plotting telemetry data and use distributions.
#' @param force Logical. If data has previously been downloaded, it will force a new
#' download and update of the OSM shapefile.
#' @param save_as_sf Logical. Should the world poolygon data be saved as
#' an \code{sf} polygon data frame? This can speed subsetting in the
#' \code{\link[crawlUtils]{cu_osm_coast}} function.
#' @references https://osmdata.openstreetmap.de/data/land-polygons.html
#' @author Josh M. London and Devin S. Johnson
#' @importFrom utils download.file unzip
#' @export
#'
# cu_download_osm <- function(force = FALSE) {
#   dir <- file.path(system.file(package="crawlUtils"), "inst", "coast")
#   if(file.exists(dir) & !force) {
#     stop("OSM data has already been downloaded. Use 'force=TRUE' to update.")
#   } else {
#     inp <- readline(prompt = "This function will download a considerable amount of coastline data.\nAre you sure you want to proceed? [y/n]: ")
#     if(tolower(inp)%in%c("n","no")) stop("Daylight OSM download crisis averted, phewww!")
#     options(timeout = max(10000, getOption("timeout")))
#     dir.create(dir, recursive=TRUE, showWarnings=FALSE)
#     download.file("https://daylight-map-distribution.s3.us-west-1.amazonaws.com/release/v1.21/coastlines-v1.21.tgz",
#                   destfile = file.path(dir,"coastlines-v1.21.tgz"),
#                   method = "auto")
#     untar(file.path(dir,"coastlines-v1.21.tgz"), exdir = dir)
#     file.remove(file.path(dir,"coastlines-v1.21.tgz"))
#     # m <- paste(readLines(file.path(dir, "land-polygons-complete-4326","README.txt")),"\n")
#     # message(m)
#   }
# }

cu_download_osm <- function(force = FALSE, save_as_sf=TRUE) {
  dir <- file.path(system.file(package="crawlUtils"), "inst", "osm")
  if(file.exists(file.path(dir,"land-polygons-complete-4326", "land_polygons.shp")) & !force) {
    stop("OSM data has already been downloaded. Use 'force=TRUE' to update.")
  } else {
    inp <- readline(prompt = "This function will download a considerable amount of coastline data.\nAre you sure you want to proceed? [y/n]: ")
    if(tolower(inp)%in%c("n","no")) stop("OSM download crisis averted, phewww!")
    options(timeout = max(1000, getOption("timeout")))
    dir.create(dir, recursive=TRUE)

    message("Downloading land polygons ...")
    download.file("https://osmdata.openstreetmap.de/download/land-polygons-complete-4326.zip",
                  destfile = file.path(dir,"land-polygons-complete-4326.zip"),
                  method = "auto")
    message("Unpacking polygons ...")
    unzip(file.path(dir,"land-polygons-complete-4326.zip"), exdir = dir)
    file.remove(file.path(dir,"land-polygons-complete-4326.zip"))

    message("Checking polygons for validity ...")
    land_path <- file.path(dir, "land-polygons-complete-4326", "land_polygons.shp")
    land <- st_read(land_path)
    chk <- st_is_valid(land)
    # all(chk)
    if(!all(chk)){
      message("Repaing nonvalid polygons ...")
      ind <- which(!chk)
      # mapview(land[ind,])
      land[ind,] <- st_make_valid(land[ind,])
    }

    if(save_as_sf){
      message("Saving {sf} polygon data ...")
      saveRDS(land, file.path(dir, "land-polygons-complete-4326", "land.rds"))
    }

    m <- paste(readLines(file.path(dir, "land-polygons-complete-4326","README.txt")),"\n")
    message(m)
  }
}

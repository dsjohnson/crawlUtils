library(crawl)
library(sf)
library(mapview); mapviewOptions(fgb=FALSE)
data(beardedSeals)

x<- st_as_sf(beardedSeals, coords=c("longitude","latitude"), crs=4326)
x <- st_transform(x, 3338)
x <- st_union(x)

download_osm()

land <- get_osm_coast(x)



% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cu_osm_coast.R
\name{cu_osm_coast}
\alias{cu_osm_coast}
\title{Get A Coastline \code{sf} Polygon Object For Plotting and Mapping}
\usage{
cu_osm_coast(x, keep = 0.2, union = TRUE)
}
\arguments{
\item{x}{An \code{sf} spatial object. The coastline will be cropped to the bounding box.}

\item{keep}{The amount of data retained after simplification with \code{\link[rmapshaper:ms_simplify]{rmapshaper::ms_simplify}}}

\item{union}{Logical. Should the retuned object be returned as a single \code{sf} \code{MULTIPOLYGON} object}
}
\description{
Uses downloaded OSM data for constructing an \code{sf} polygon coastline data object.
Prior to using this function you must run \code{crawlUtils::cu_download_osm()}.
See \code{\link[crawlUtils:cu_download_osm]{crawlUtils::cu_download_osm}}.
}
\author{
Josh M. London and Devin S. Johnson
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sf_funcs.R
\name{st_to_360}
\alias{st_to_360}
\title{Convert longitude/latitude coordinates from -180/180 to 0/360}
\usage{
st_to_360(x)
}
\arguments{
\item{x}{An sf data frame with EPSG=4326.}
}
\description{
Converts sf data with EPSG = 4326 from -180/180 specification to
0/360 for plotting with the mapview package etc.
}
\author{
Josh London
}

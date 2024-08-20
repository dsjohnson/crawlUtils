#' @title Movement bout detection
#' @description Creates a data table that indicates the times of different
#' bouts of movement. This method uses changes in the overall dispersion rate of the animal
#' from the 'base' time to detect changes in overall movement from small scale local movement
#' to large scale migration.
#' @param bout_tbl A bout table for an individual animal created by \code{\link[crawlUtils]{cu_bout_det}}
#' @param data An \code{sf} data set of locations with times noted by the 'datetime' column
#' @param grid_res The temporal resolution at which migrations are detected. e.g., "day" (default) implies migration
#' start and end is detected on a daily resolution.
#' @param base The location at which dispersion is measured. Can be one of \code{"first"} (first location),
#' \code{"last"} (final location), or some other \code{sf::sfc} point location.
#' @export
#' @author Devin S. Johnson
#' @import units lubridate dplyr sf
#' @importFrom stats coef dist predict vcov coefficients lm
#' @importFrom tidyr nest
#'
cu_disp_rate <- function(bout_tbl, data, grid_res="day", base="first"){
  travel <- bout <- datetime <- avg_disp_rate <- i <- NULL
  bout_tbl <- bout_tbl %>% select(-avg_disp_rate)
  if(base=="first"){
    base <- data[1,] %>% st_geometry()
  } else if(base=="last"){
    base <- data[nrow(data),] %>% st_geometry()
  } else if(!inherits(base,c("sf","sfc"))){
    stop("Unrecognized 'base' argument!")
  }
  ddd <- data.frame(
    datetime = data$datetime,
    dist = st_distance(data, base) %>% units::set_units("km") %>% as.vector()
  )
  ddd$time <- with(ddd,
                   as.numeric(datetime)/as.numeric(duration(1, grid_res))
  )
  ddd <- cu_add_bouts(ddd, bout_tbl) %>% select(datetime, bout, dist)
  ddd <- ddd %>% group_by(bout) %>% nest()
  ddd$avg_disp_rate <- foreach(i=1:nrow(ddd), .combine = c)%do%{
    rate <-  as.numeric(max(ddd$data[[i]]$dist) - min(ddd$data[[i]]$dist))/as.numeric(diff(range(ddd$data[[i]]$datetime)))
    units::set_units(rate,  paste0("km/",grid_res), mode='standard')
  }
  ddd$total_disp <- foreach(i=1:nrow(ddd), .combine = c)%do%{
    max(ddd$data[[i]]$dist) - min(ddd$data[[i]]$dist)
  }
  ddd <- ddd %>% select(-data)

  ddd <- full_join(ddd, bout_tbl, by="bout")

  attr(ddd, "base") <- base

  return(ddd)
}

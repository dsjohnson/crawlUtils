#' @title Label telemetry bout segments
#' @description Label segments of telemetry data such that gaps between observed
#' locations are less than or equal to the interval specified.
#' @param x A data set of telemetry locations and times.
#' @param gap specified maximum time interval for which a new bout is started. Defaults
#' to \code{gap=7} days.
#' @param time_unit Unit of time of the gap specification. Defaults to \code{"day"}.
#' @param ... Ignored arguments
#' @import lubridate
#' @author Devin S. Johnson
#' @export
#'
cu_add_gaps <- function(x, gap=7, time_unit="days", ...){
  . <- NULL #handle 'no visible binding...'
  time <- x$datetime
  dt <- diff(time) %>% `units<-`(time_unit)
  time_diff <- c(0, dt)
  x$bout_id <- (time_diff >= gap) %>% {cumsum(.)+1}
  return(x)
}

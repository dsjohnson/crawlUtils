#' @title Label telemetry bout segments
#' @description Label segments of telemetry data such that gaps between observed
#' locations are less than or equal to the interval specified.
#' @param x column of POSIX times for which locations are observed.
#' @param gap specified maximum time interval for which a new bout is started. Defaults
#' to \code{gap=7} days.
#' @param time_unit Unit of time of the gap specification. Defaults to \code{"day"}.
#' @import lubridate
#' @author Devin S. Johnson
#' @export
#'
cu_get_bouts <- function(x, gap=7, time_unit="days"){
  . <- NULL #handle 'no visible binding...'
  dt <- diff(x) %>% `units<-`(time_unit)
  time_diff <- c(0, dt)
  bout_id <- (time_diff >= gap) %>% {cumsum(.)+1}
  return(bout_id)
}

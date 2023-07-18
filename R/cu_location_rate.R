#' @title Location rate statistics
#' @description Calculate location rate statistics such as mean location rate per
#' day or maximum number of locations per day for a telemetry
#' @param x data set containing time of locations.
#' @param time_name Character name of the POSIX time column when locations are
#' observed.
#' @param time_unit Time unit of the location summary. Defaults to \code{"day"}.
#' @param stat Function used to summarize location times. Defaults to \code{mean}.
#' @param ... Additional arguments passed to \code{stat} function.
#' @author Devin S. Johnson
#' @import lubridate
#' @export
#'
cu_location_rate <- function(x, time_name, time_unit="day", stat=mean, ...){
  stop("This function is currently broken!")
  out <- table(lubridate::round_date(x[,time_name], time_unit))
  return(stat(out,...))
}

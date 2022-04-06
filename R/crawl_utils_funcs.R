

#' @title Label telemetry bout segments
#' @description Label segments of telemetry data such that gaps between observed
#' locations are less than or equal to the interval specified.
#' @param x column of POSIX times for which locations are observed.
#' @param gap specified maximum time interval for which a new bout is started. Defaults
#' to \code{gap=7} days.
#' @param time_unit Unit of time of the gap specification. Defaults to \code{"day"}.
#' @import lubridate units
#' @author Devin S. Johnson
#' @export
#'
cu_get_bouts <- function(x, gap=7, time_unit="day"){
  . <- NULL #handle 'no visible binding...'
  dt <- diff(x) %>% `units<-`(time_unit)
  time_diff <- c(0, dt)
  bout_id <- (time_diff >= gap) %>% {cumsum(.)+1}
  return(bout_id)
}


#' @title Join crawl prediction or simulation output with a table based on
#' a time interval
#' @description Takes a data set with a POSIX time column named 'datetime' and
#' another data set with \code{start} and \code{end} columns representing time intervals and
#' merges the two depending whether or not the 'datetime' column is within the interval  of the second.
#' @param x A data frame with a column labeled \code{datetime}
#' @param int_tbl A data frame with 'start' and 'end' columns that form non-overlapping intervals as well as at least
#' one other column with interval level data.
#' @author Devin S. Johnson
#' @importFrom fuzzyjoin fuzzy_left_join
#' @export
#'
cu_join_interval_tbl <- function(x, int_tbl){
  x <- fuzzyjoin::fuzzy_left_join(x,int_tbl,
                                  by=c(datetime="start",datetime="end"),
                                  match_fun = list(`>=`, `<=`)
  )
}


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

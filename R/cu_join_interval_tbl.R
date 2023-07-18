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

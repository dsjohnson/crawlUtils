#' @title Create a Table of Start and End Times for each Bout in a telemetry data set
#' @description Summarizes a telemetry dataframe with a bout column into
#' a data frame with start and end times for each bout. For use with
#' \code{\link[crawlUtils:cu_join_interval_tbl]{cu_join_interval_tbl}}.
#' @param x A telelmtry data frame containing a bout column
#' @param bout Character. The name of the bout column.
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr group_by summarize
#' @author Devin S. Johnson
#' @export
cu_bout_summary <- function(x,bout){
  datetime <- NULL
  x <- x %>% group_by(.data[[bout]]) %>% st_drop_geometry() %>%
    summarize(start=min(datetime), end=max(datetime))
}

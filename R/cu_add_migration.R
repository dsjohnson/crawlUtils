#' @title Add migration detection results to original location data
#' @description Add bout data from \code{migr_evt}
#' See \code{\link[crawlUtils]{cu_bout_det}}
#' @param data Original data used by \code{cu_migration_det}.
#' @param bout_tbl Results table produced by \code{cu_migration_det}.
#' @export
#' @author Devin S. Johnson
#' @import sf fuzzyjoin dplyr ggplot2

cu_add_bouts <- function(data, bout_tbl){
  travel <- bout <- datetime <- start <- end <- avg_disp_rate <- NULL
  if("travel" %in% colnames(data)) data <- select(data, -travel)
  if("bout" %in% colnames(data)) data <- select(data, -bout)
  bout_end <- bout_tbl$start[-1]
  bout_tbl <- head(bout_tbl,-1)
  bout_tbl$end <- bout_end
  data <- fuzzyjoin::fuzzy_left_join(data,bout_tbl,
                                     by=c(datetime="start",datetime="end"),
                                     match_fun = list(`>=`, `<`))
  data <- arrange(data, datetime)
  return(data)
}

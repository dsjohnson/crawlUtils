#' @title Add migration detection results to original location data
#' @description Add \code{migr_evt}
#' See \code{\link[crawlUtils]{cu_migration_det}}
#' @param data Original data used by \code{cu_migration_det}.
#' @param migr_tbl Results table produced by \code{cu_migration_det}.
#' @export
#' @author Devin S. Johnson
#' @import sf fuzzyjoin dplyr ggplot2

cu_add_migration <- function(data, migr_tbl){
  migr_evt <- phase <- datetime <- start <- end <- avg_disp_rate <- NULL
  if("migr_evt" %in% colnames(data)) data <- select(data, -migr_evt)
  if("phase" %in% colnames(data)) data <- select(data, -phase)
  data <- fuzzyjoin::fuzzy_left_join(data,migr_tbl,
                                     by=c(datetime="start",datetime="end"),
                                     match_fun = list(`>=`, `<`))
  starts <- migr_tbl %>% select(start, migr_evt, phase) %>%
    rename(datetime=start) %>% slice(-1)
  data <- full_join(data, starts, by = c("datetime", "migr_evt", "phase"))
  data <- data %>% select(-start, -end, -avg_disp_rate) %>%
    arrange(datetime)
  return(data)
}

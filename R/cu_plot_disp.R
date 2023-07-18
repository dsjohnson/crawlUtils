#' @title Plot migration detection results
#' @description PLot the animal dispersion from the base location over time. Points
#' are colored to reflect estimated migration and non-migration phases.
#' See \code{\link[crawlUtils]{cu_migration_det}}
#' @param data Original data used by \code{cu_migration_det}.
#' @param migr_tbl Results table produced by \code{cu_migration_det}.
#' @export
#' @author Devin S. Johnson
#' @import sf fuzzyjoin dplyr ggplot2
cu_plot_disp <- function(data, migr_tbl){
  datetime <- dist <- migr_evt <- NULL
  base <- attr(migr_tbl, "base")
  ddd <- data.frame(
    datetime = data$datetime,
    dist = st_distance(data, base) %>% units::set_units("km") %>% as.vector()
  )
  ddd <- fuzzyjoin::fuzzy_left_join(ddd,migr_tbl,
                                    by=c(datetime="start",datetime="end"),
                                    match_fun = list(`>=`, `<=`))
  plt <- ggplot() +
    geom_point(aes(x=datetime, y=dist), alpha=1, color="slategray2",
               data=ddd %>% filter(migr_evt=="0")) +
    geom_point(aes(x=datetime, y=dist), alpha=1, color="darkred",
               data=ddd %>% filter(migr_evt=="1")) +
    theme_classic() +
    xlab("Date") + ylab("Dispersal (km)")
  return(plt)
}

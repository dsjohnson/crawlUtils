#' @title Plot migration detection results
#' @description PLot the animal dispersion from the base location over time. Points
#' are colored to reflect estimated migration and non-migration phases.
#' See \code{\link[crawlUtils]{cu_bout_det}}
#' @param data Original data used by \code{cu_bout_det}.
#' @param bout_tbl Results table produced by \code{cu_bout_det}.
#' @export
#' @author Devin S. Johnson
#' @import sf fuzzyjoin dplyr ggplot2
cu_plot_disp <- function(data, bout_tbl){
  datetime <- dist <- travel <- NULL
  base <- attr(bout_tbl, "base")
  ddd <- data.frame(
    datetime = data$datetime,
    dist = st_distance(data, base) %>% units::set_units("km") %>% as.vector()
  )
  bout_end <- bout_tbl$start[-1]
  bout_tbl <- head(bout_tbl,-1)
  bout_tbl$end <- bout_end
  ddd <- fuzzyjoin::fuzzy_left_join(ddd,bout_tbl,
                                    by=c(datetime="start",datetime="end"),
                                    match_fun = list(`>=`, `<=`))
  plt <- ggplot() +
    geom_point(aes(x=datetime, y=dist), alpha=1, color="slategray2",
               data=ddd %>% filter(travel=="0")) +
    geom_point(aes(x=datetime, y=dist), alpha=1, color="darkred",
               data=ddd %>% filter(travel!="0")) +
    theme_classic() +
    xlab("Date") + ylab("Dispersal (km)")
  return(plt)
}

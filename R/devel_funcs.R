#' @title Convert sf output from crawl to ctmm object
#' @description Takes \code{crawl::crwPredict} or \code{crawl::crwPostIS}
#' objects converted to \code{sf} objects by \code{\link[crawl:crw_as_sf]{crawl::crw_as_sf} and converts them to data
#' objects that can be used within the \code{ctmm} package.
#' @param data \code{sf} data object from \code{\link[crawl:crw_as_sf]{crawl::crw_as_sf}
#' @param silent Logical. Should messages from \code{ctmm} be printed. Defaults to
#' \code{silent = FALSE}.
#' @author Devin S. Johnson
#' @import sf dplyr magrittr
#' @importFrom
#' @export
#'
cu_crw2ctmm <- function(data, silent=FALSE){
  warning("'crw2ctmm()' is highly experimental, use is not recommended at this time. Use with caution!")
  if (! requireNamespace("ctmm", quietly = TRUE)) stop("Please install {ctmm}: install.packages('ctmm')")
  datetime <- NULL
  tmp <- data  %>% sf::st_transform(4326) %>% sf::st_coordinates() %>%
    cbind(data,.) %>% dplyr::select(datetime, X, Y, rep) %>% sf::st_drop_geometry() %>%
    mutate(timestamp=as.character(datetime)) %>%
    rename(
      `location-long` = X,
      `location-lat` = Y,
      `tag-local-identifier` = rep
    )
  if(silent){
    out <- suppressMessages(
      ctmm::as.telemetry(
        tmp,
        projection=sf::st_crs(data)$proj4string)
    )
  } else{
    out <- ctmm::as.telemetry(
      tmp,
      projection=sf::st_crs(data)$proj4string
    )
  }
  return(out)
}

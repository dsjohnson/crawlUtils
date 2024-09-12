#' Calculate an empirical utilization distribution
#'
#' @param pts an sf point object
#' @param grid a  `SpatRaster` on which to evaluate the UD.
#' @param average should an average UD be created from 1 or more pts layers
#'
#' @return grid with an additional column, _npts_ or _mean_pts_
#' @export
#'
cu_empirical_ud <- function(pts, grid, average = TRUE) {
  cell <- npts <- NULL
  npts_in_poly <- function(x, y) {
    y$npts = lengths(sf::st_intersects(y, x))
    return(y)
  }

  if(inherits(pts,'sf')) {
    pts <- list(pts)
  }
  res <- purrr::map(pts, ~npts_in_poly( x = ., y = grid))
  if(average) {
    res <- res %>%
      dplyr::bind_rows() %>%
      sf::st_drop_geometry() %>%
      dplyr::group_by(cell) %>%
      dplyr::summarise(mean_pts = mean(npts))
    res <- grid %>% dplyr::left_join(res, by="cell")
    return(res)
  }
  return(purrr::pluck(res,1))
}

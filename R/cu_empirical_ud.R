#' Calculate an empirical utilization distribution
#'
#' @param pts an sf point object
#' @param grid a grid created from [cu_ud_grid()]
#' @param average should an average UD be created from 1 or more pts layers
#'
#' @return grid with an additional column, mean_pts
#' @export
#'
cu_empirical_ud <- function(pts, grid, average = TRUE) {
  
  npts_in_poly <- function(x, y) {
    res <- dplyr::mutate(y, npts = lengths(sf::st_intersects(y, x)))
  }
  
  if(inherits(pts,'sf')) {
    pts <- list(pts)
  }
  res <- purrr::map_df(pts, npts_in_poly, y = grid)
  
  if(average) {
    res <- res %>% 
      dplyr::bind_rows() %>% 
      sf::st_drop_geometry() %>% 
      dplyr::group_by(cell) %>% 
      dplyr::summarise(mean_pts = mean(npts))
    res <- grid %>% dplyr::left_join(res, by="cell")
  }
  return(res)
}

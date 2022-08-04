#' Calculate an empirical utilization distribution
#'
#' @param grid 
#' @param pts 
#' @param average 
#'
#' @return
#' @export
#'
#' @examples
cu_empirical_ud <- function(grid, pts, average = TRUE) {
  npts_in_poly <- function(grid, pts) {
    res <- grid %>% 
      dplyr::mutate(npts = lengths(sf::st_intersects(grid, pts)))
  }
  res <- purrr::map2(list(grid), pts, npts_in_poly)
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

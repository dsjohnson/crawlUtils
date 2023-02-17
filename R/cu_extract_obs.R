#' @title Extract predicted locatiins closest in time to observed locations.
#' @param x A `crwPredict` or `crwIS` object. Time column must be named `datetime`.
#' @param obs Observed locations. Time column must be named `datetime`.
#' @description This function is intended for use after a predicted track or posterior simulation
#' is routed around barriers using the `{pathroutr}` package.
#' @author Devin S. Johnson
#' @export
cu_extract_obs <- function(x, obs){
  a <- obs$datetime
  b <- x$datetime
  ind <- apply(outer(a, b, \(x,y) abs(x-y)), 1, which.min)
  return(x[ind,])
}

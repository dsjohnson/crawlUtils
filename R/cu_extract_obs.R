#' @title Extract predicted locatiins closest in time to observed locations.
#' @param x A `crwPredict` or `crwIS` object. Time column must be named `datetime`.
#' @param obs Observed locations. Time column must be named `datetime`.
#' @description This function is intended for use after a predicted track or posterior simulation
#' is routed around barriers using the `{pathroutr}` package.
#' @author Devin S. Johnson
#' @export
cu_extract_obst <- function(x, obs){
  a <- obs$datetime
  x_type <- attr(x, "crw_type")
  if(x_type %in% c("crwIS_list", "crwIS_sf_list")){
    x2_type <- attr(x[[1]], "crw_type")
    b <- x[[1]]$datetime
  } else if(x_type %in% c("crwPredict","crwPredict_sf", "crwIS","crwIS_sf")){
    b <- x$datetime
  }else{
    stop("x has unknown attr(x, 'crw_type')")
  }
  ind <- apply(outer(a, b, \(x,y) abs(x-y)), 1, which.min)
  if(x_type %in% c("crwPredict","crwPredict_sf", "crwIS","crwIS_sf")){
    x <- x[ind,]
    attr(x, "crw_type") <- x_type
  } else{
    for(i in 1:length(x)){
      x[[i]] <- x[[i]][ind,]
      attr(x[[i]], "crw_type") <- x2_type
    }
    attr(x, "crw_type") <- x_type
  }
return(x)
}

#' @title Compute velocity based kernel density estimate bandwidth
#' @param pts A prediction or simulation data set produced by \code{\link[crawlUtils]{cu_crw_predict}} or
#' \code{\link[crawlUtils]{cu_crw_sample}}
#' @param ess An effective sample size value
#' @param vel_quant A quantile value for the observed average velocities, e.g., 0.5 would imply use of the median
#' observed average velocity between time points
#' @param vel_fix A fixed value for the average velocity used
#' @param time_scale The time scale used, e.g., movement per \code{"hours"}.
#' @export
#' @importFrom lubridate duration
#' @importFrom stats quantile
cu_vel_B <- function(pts, ess=nrow(pts), vel_quant=0.5, vel_fix=NULL, time_scale="hours"){
  if(!inherits(pts, "sf")) stop("Input 'pts' needs to be of class 'sf'!")
  time <- pts$datetime
  dt <- as.numeric(diff(time))/as.numeric(lubridate::duration(1,time_scale))
  if(!is.numeric(vel_fix)){
    xy <- st_coordinates(pts)
    step <- sqrt(diff(xy[,1])^2 + diff(xy[,2])^2)
    avg_vel <- step/dt
    vel_fix <- quantile(avg_vel, prob=vel_quant, na.rm=TRUE)
  }
  rt <- diff(range(time))
  units(rt) <- time_scale
  tpo <- as.numeric(rt/ess)
  B <- diag(rep((vel_fix*tpo/2)^2,2))
  return(B)
}

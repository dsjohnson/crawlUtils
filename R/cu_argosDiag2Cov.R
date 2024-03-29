#' @title Transform Argos diagnostic data to covariance matrix form
#'
#' @description Using this function the user can transform the Argos diagnostic data for location
#' error into a form usable as a covariance matrix to approximate the location error with a
#' bivariate Gaussian distribution. The resulting data.frame should be attached back to the data
#' with \code{cbind} to use with the \code{crwMLE} function.
#' @param Major A vector containing the major axis information for each observation (na values are ok)
#' @param Minor  A vector containing the minor axis information for each observation (na values are ok)
#' @param Orientation A vector containing the angle orientation of the Major axis from North (na values are ok)
#' @return A \code{data.frame} with the following columns
#' \item{ln.sd.x}{The log standard deviation of the location error in the x coordinate}
#' \item{ln.sd.y}{The log standard deviation of the location error in the x coordinate}
#' \item{rho}{The correlation of the bivariate location error ellipse}
#' @author Devin S. Johnson
#' @export

cu_argosDiag2Cov = function(Major, Minor, Orientation){
  a=Major
  b=Minor
  if(any(b<=.Machine$double.eps, na.rm=TRUE)) stop("There are very small (or 0) values for the minor ellipse lengths! These may need to be removed.")
  theta=Orientation
  if(any(theta < 0 | theta > 180, na.rm = TRUE)) stop("Argos diagnostic data orientation outside of [0,180]!")
  if(any(a < 0, na.rm = TRUE)) stop("Argos diagnostic data major axis < 0!")
  if(any(b < 0, na.rm = TRUE)) stop("Argos diagnostic data minor axis < 0!")
  theta = pi*(theta/180)
  k=sqrt(2)
  v1 = (a/k)^2*sin(theta)^2 + (b/k)^2*cos(theta)^2
  v2 = (a/k)^2*cos(theta)^2 + (b/k)^2*sin(theta)^2
  c12 = ((a^2 - b^2)/k^2)*cos(theta)*sin(theta)
  rho = c12/(sqrt(v1)*sqrt(v2))
  check = (v1*v2-c12^2) > 0
  if(any(rho > 1 | rho < -1, na.rm=TRUE)) stop("Faulty Argos error correlation calculated from 'argosDiag2Cov' function")
  return(data.frame(ln.sd.x=log(sqrt(v1)), ln.sd.y=log(sqrt(v2)), error.corr=rho, diag.check=check))
}

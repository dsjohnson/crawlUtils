#' @title Create covariance function for a fitted CRW model object
#' @description A function is created to evaluate the covariance function of the fitted
#' CRW movement model
#' @param x A crwFit object created by a call to \code{\link[crawl]{crwMLE}}
#' @details The function returns a function to evaluate the covariance of the fitted Integrated
#' Ornstein-Ulenbeck movement model. The returned function has 3 arguments: (1) \code{t1}
#' and (2) \code{t2} both vectors of times to evaluate the covariance function of the fitted
#' IOU model, and (3) \code{E} (defaults to \code{E=0}). Which is the "zero" time of the process.
#' Typically \code{E} will be the time of the first observation.
#' @references Taylor, J. M., Cumberland, W. G., & Sy, J. P. (1994).
#' A stochastic model for analysis of longitudinal AIDS data. Journal of the
#' American Statistical Association, 89(427), 727-736.
#'
#' Johnson, D. S., London, J. M., Lea, M. A., & Durban, J. W. (2008).
#' Continuous‐time correlated random walk model for animal telemetry data.
#' Ecology, 89(5), 1208-1215.
#'@export
#'@importFrom utils tail
#'@importFrom stats var
#'@author Devin S. Johnson
#'
cu_crw_covfun <- function(x){
  if((x$mov.model != ~1) | (x$random.drift==TRUE) | (!is.null(x$activity))) stop("Sorry, currently this function only works with the base CRW model")
  par <- tail(x$par,2)
  b <- exp(par[2])
  sig2 <- exp(2*par[1])
  Pvec <- diag(var(x$data[,x$coord]))
  foo <- function(t1, t2, E=min(x$data$TimeNum)){
    if(E>min(c(t1,t2))) stop("E must be < min(c(t1,t2))")
    s <- pmin(t1-E, t2-E)
    t <- pmax(t1-E, t2-E)
    sig2*(s - (1/(2*b))*(1+exp(-b*(t-s))) + (1/(2*b))*exp(-b*t) + (1/(2*b))*exp(-b*s) )
  }
  attr(foo, "b") <- b
  attr(foo, "sig2") <- sig2
  return(foo)
}

#' @title Calculate correlation matrix for a set of times from a CRW covariance function
#' @description Using a correlation function created by \code{\link{cu_crw_covfun}}
#' from a fitted CRW model a covariance (correlation) matrix is created for observations
#' at the user provided times.
#' @param x Either a \code{crwFit} object from a call to \code{\link[crawl]{crwMLE}}, or a
#' vector of times.
#' @param corr Should the function return a correlation or covariance matrix? Defaults
#' to \code{corr = TRUE}.
#' @param cf Covariance function created from a \code{crwFit} object with a call to
#' \code{\link{cu_crw_covfun}}.
#' @param E The 'zero' time used for the covariance function. Defaults to \code{E = 0}.
#' (See \code{\link{cu_crw_covfun}}).
#' @details If \code{x} is a \code{crwFit} object, then the \code{cf} and \code{E}
#' arguments are ignored. The resulting matrix is the covariance (correlation) matrix
#' for the observed location times conditioned on the first location. Therefore,
#' for n observed locations an (n-1) by (n-1) matrix will result. This is most
#' useful for the effective sample size computation for a kernel density estimate.
#' @author Devin S. Johnson
#' @importFrom stats cov2cor
#' @export
cu_crw_covmat <- function(x, corr=TRUE, cf, E=0){
  if(inherits(x, "crwFit")){
    foo <- cu_crw_covfun(x)
    E <- min(x$data$TimeNum)
    times <-
      S <- outer(x$data$TimeNum[-1],x$data$TimeNum[-1], foo, E=E)
  } else {
    S <- outer(x,x,FUN=cf,E=E)
  }
  if(corr){
    S <- cov2cor(S)
  }
  return(S)
}

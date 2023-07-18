#' @title Batch Fitting CTCRW Models for Argos (and FastGPS) Data
#' @description A basic CTCRW model is fitted to a list of data sets where each
#' element in the list represents the data for a single individual or deployment.
#' @param data A data set containg the telemetry data for a single individual.
#' Will also accept a single \code{sf} data frame as well.
#' @param move_phase An optional character value indicating a factor variable column in the data that
#' designates different movement phases.
#' @param bm Fit a Brownian Motion model rather than in integrated OU model. Defaults to \code{bm = FALSE}.
#' @param use_prior Logical. Should a sensible mixture normal prior be use for the log beta and
#' log error scale parameters to impose a soft constraint for better numerical optimization.
#' Default is \code{TRUE}.
#' @param crw_control A named list passed to \code{\link[crawl]{crwMLE}} for optimization. Alternatives
#' for the default values of \code{initialSANN}, \code{attempts}, \code{control}, \code{theta}, \code{fixPar}, \code{constr}, and \code{prior}
#' can be specified here. See \code{\link[crawl]{crwMLE}} for a description of these arguments. WARNING!!! No
#' checks are made for validity of the user override. So know what you are doing.
#' @param fit Logical. CTCRW parameters are estimated if \code{fit=TRUE} (default), else
#' the results of \code{\link[crawl]{displayPar}}.
#' @param skip_check See \code{\link[crawl]{crwMLE}} v2.3.0. Currently ignored.
#' @param ... Additional arguments passed to the \code{\link[foreach]{foreach}} function, e.g.,
#' for error handling in the loop.
#' @import dplyr crawl sf foreach
#' @importFrom stats as.formula dexp model.frame model.matrix na.pass qchisq
#' @importFrom utils head
#' @importFrom stats dnorm
#' @export
#'
cu_crw_argos <- function(data, move_phase=NULL, bm=FALSE, use_prior=TRUE,
                         crw_control=NULL, fit = TRUE, skip_check=FALSE,...){
  datetime <- type <- const <- NULL #handle 'no visible binding...'

  if(!all(c('ln.sd.x','ln.sd.y','error.corr')%in%colnames(data))){
    stop("It appears that argos columns have not been added via 'cu_add_argos_cols()'")
  }

  # progressr::handlers(global = TRUE)
  # if(length(data_list)>1) p <- progressr::progressor(length(data_list))
  # Define mix normal priors:
  if(use_prior){
    pb <- function(x){log(mean(dnorm(x,seq(-8.5,1.5,0.5),0.5)))}
    prior_b <- function(x){sum(sapply(x,pb))}
    plq <- function(x){log(mean(dnorm(x,seq(0,0.5,0.1),0.1)))}
    prior_lq <- function(x){sum(sapply(x,plq))}
  }

  if (inherits(data,"sf")) {
    proj_unit <- sf::st_crs(data, parameters = TRUE)$units_gdal
    if (proj_unit == "metre") {
      units <- 1
    } else if (proj_unit == "kilometre") {
      units <- 1/1000
    } else{
      if(units %in% c('meter','metre')) units <- 1
      if(units %in% c('kilometer','kilometre')) units <- 1/1000
    }
  } else {
    units <- 1
  }

  data <- data |> dplyr::arrange(datetime)

  if(!is.null(move_phase)){
    mov.model <- as.formula(paste0("~0+",move_phase))
    mov.mf <-
      model.matrix(mov.model, model.frame(mov.model, data, na.action = na.pass))
    if (any(is.na(mov.mf)))
      stop("Missing values are not allowed in move_phase!")
    n.mov <- ncol(mov.mf)
  }else{
    n.mov <- 1
    mov.model <- ~1
  }

  # movement parameters quantities:
  if(!bm){
    mov.theta <- c(rep(8+log(units),n.mov),rep(log(-log(0.1)),n.mov))
    mov.fix <- rep(NA, 2*n.mov)
    if(use_prior){
      mov.prior <- function(par){prior_b(tail(par,n.mov))}
    } else{
      mov.prior <- function(par){return(0)}
    }

  } else if(bm){
    mov.theta <- c(rep(8 + log(units),n.mov))
    mov.fix <- c(rep(NA, n.mov), rep(3, n.mov))
    mov.prior <- function(par){return(0)}
  } else{
    stop("The 'bm' argument should be TRUE or FALSE")
  }

  err.model <- list(x=~ln.sd.x, y=~ln.sd.y, rho= ~error.corr)
  err.fix <- c(NA,1,NA,1)
  err.theta <- c(0,0)
  #err.prior <- function(par){return(0)} #function(par){prior_lq(par[1:2])}
  if(use_prior){
    err.prior <- function(par){prior_lq(par[1:2])}
  } else{
    err.prior <- function(par){return(0)}
  }

  # Fit ctcrw model
  if(is.null(crw_control$initialSANN)){
    initialSANN <- list(maxit=1500, temp=10)
  }else{
    initialSANN <- crw_control$initialSANN
  }
  if(is.null(crw_control$attempts)){
    attempts <- 20
  }else{
    attempts <- crw_control$attempts
  }
  if(is.null(crw_control$control)){
    control <- list(maxit=10000)
  }else{
    control <- crw_control$control
  }
  if(is.null(crw_control$theta)){
    theta <- c(err.theta, mov.theta)
  }else{
    theta <- crw_control$theta
  }
  if(is.null(crw_control$fixPar)){
    fixPar <- c(err.fix, mov.fix)
  } else{
    fixPar <- crw_control$fixPar
  }
  if(is.null(crw_control$prior)){
    prior <- function(par){err.prior(par) + mov.prior(par)}
  } else{
    prior <- crw_control$prior
  }
  if(is.null(crw_control$method)){
    method <- "Nelder-Mead"
  } else{
    method <- crw_control$method
  }


  if(fit){
    if(is.null(crw_control$constr)){
      constr <- list(lower = -Inf, upper = Inf)
      suppressMessages(
        out <- crawl::crwMLE(
          mov.model = mov.model, err.model = err.model, data = data, Time.name="datetime",
          fixPar = fixPar, theta = theta,
          control = control, initialSANN = initialSANN,
          prior=prior,
          attempts=attempts, method = method)
      )
    } else{
      suppressMessages(
        out <- crawl::crwMLE(
          mov.model = mov.model, err.model = err.model, data = data, Time.name="datetime",
          fixPar = fixPar, theta = theta,
          control = control, constr=crw_control$constr, initialSANN = initialSANN,
          prior=prior,
          attempts=attempts, method = "L-BFGS-B")
      )
    }
  } else{
    suppressMessages(
      out <- crawl::displayPar(
        mov.model = mov.model, err.model = err.model, data = data, Time.name="datetime",
        fixPar = fixPar, theta = theta,
        control = control, constr=constr, initialSANN = initialSANN,
        prior=prior,
        attempts=attempts)
    )
  }
  attr(out, "crw_type") <- "crwFit_argos"
  out
  return(out)
}

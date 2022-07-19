#' @title Add columns for modeling FastGPS and ARGOS error structure in the same
#' telemetry deployment
#' @description Columns are added to the telemetry data set so that multiple data types
#' can be used simultaneously within the same animal: FastGPS, Argos least-squares, and
#' Argos Kalman filter.
#' @param x Data frame containing location telemetry data and Argos quality
#' information. See 'Details' for a description of the necessary data column names.
#' @details To use this function the data set must contain the following columns with
#' exact names: (1) \code{"type"}, indicate the type of location,
#' (2) \code{"quality"} which indicates the Argos quality of the location, (3) The
#' Argos KF diagnostics: \code{"error_semi_major_axis"}, \code{"error_semi_minor_axis"}, and
#' \code{"error_ellipse_orientation"}. If there are no Argos KF locations, the associated columns
#' need not be present. Values of \code{type} are 'FastGPS', 'Argos', or 'known'. For 'FastGPS' an error radius
#' of 100m is assumed. For 'known' an error radius of 20m is assumed. Values of
#' \code{quality} must be '3','2','1','0','A',or 'B' for \code{type=='Argos'} locations and '4'-'10' for \code{type=='FastGPS'}.
#' For other types it is ignored.
#' @author Devin S. Johnson
#' @import dplyr crawl
#' @export
#'
cu_add_argos_cols <- function(x){
  error.corr <- quality <- ln.sd.x <- ln.sd.y <- error_area <- NULL #handle 'no visible binding...'
  col_nms <- colnames(x)
  if('error_semi_major_axis' %in% col_nms &
     'error_semi_minor_axis' %in% col_nms &
     'error_ellipse_orientation' %in% col_nms){
    x <- cbind(x,
               crawl::argosDiag2Cov(
                 x$error_semi_major_axis,
                 x$error_semi_minor_axis,
                 x$error_ellipse_orientation)
    )
  } else{
    x$ln.sd.x <- NA_real_; x$ln.sd.y <- NA_real_; x$error.corr <- 0
  }
  kf_ind <- !(is.na(x$ln.sd.x) | is.na(x$ln.sd.y) | is.na(x$error.corr))
  x <- x %>%
    dplyr::mutate(
      type = case_when(
        type %in% c("Argos","argos","KF") & kf_ind ~ "Argos_kf",
        type %in% c("Argos","argos","LS") & !kf_ind ~ "Argos_ls",
        type %in% c("FastGPS","GPS","G","gps") ~ "FastGPS",
        TRUE ~ type
      ),
      ln.sd.x = dplyr::case_when(
        type == "known" ~ log(10),
        type=="FastGPS" & quality=="4" ~ log(1163/2),
        type=="FastGPS" & quality=="5" ~ log(169/2),
        type=="FastGPS" & quality=="6" ~ log(71/2),
        type=="FastGPS" & quality=="7" ~ log(43/2),
        type=="FastGPS" & quality=="8" ~ log(34/2),
        type=="FastGPS" & quality=="9" ~ log(28/2),
        type=="FastGPS" & quality=="10" ~ log(24/2),
        type=="FastGPS" & quality=="11" ~ log(19/2),
        type=="Argos_ls" & quality=="3" ~ log(250),
        type=="Argos_ls" & quality=="2" ~ log(500),
        type=="Argos_ls" & quality %in% c("1","0","A","B") ~ log(1500),
        TRUE ~ ln.sd.x
      ),
      ln.sd.y = dplyr::case_when(
        type == "known" ~ log(10),
        type=="FastGPS" & quality=="4" ~ log(1163/2),
        type=="FastGPS" & quality=="5" ~ log(169/2),
        type=="FastGPS" & quality=="6" ~ log(71/2),
        type=="FastGPS" & quality=="7" ~ log(43/2),
        type=="FastGPS" & quality=="8" ~ log(34/2),
        type=="FastGPS" & quality=="9" ~ log(28/2),
        type=="FastGPS" & quality=="10" ~ log(24/2),
        type=="FastGPS" & quality=="11" ~ log(19/2),
        type=="Argos_ls" & quality=="3" ~ log(250),
        type=="Argos_ls" & quality=="2" ~ log(500),
        type=="Argos_ls" & quality %in% c("1","0","A","B") ~ log(1500),
        TRUE ~ ln.sd.y
      ),
      error.corr = ifelse(is.na(.data$error.corr), 0, .data$error.corr),
      gq5 = ifelse(quality %in% c("4","5"), 1, 0),
      gq4 = ifelse(quality=="4", 1, 0),
      aq0 = ifelse(quality %in% c("0","A","B"), 1, 0),
      aqA = ifelse(quality %in% c("A","B"), 1, 0),
      aqB = ifelse(quality=="B", 1, 0),
      error_area = pi*qchisq(0.95, 2)*exp(0.5*(ln.sd.x+ln.sd.y))*sqrt(1-error.corr^2),
      error_area = ifelse(quality=="A", 1.1*error_area, error_area),
      error_area = ifelse(quality=="B", 1.2*error_area, error_area)
    )
  return(x)
}



#' @title Batch Fitting CTCRW Models for Argos (and FastGPS) Data
#' @description A basic CTCRW model is fitted to a list of data sets where each
#' element in the list represents the data for a single individual or deployment.
#' @param data_list A list of data sets. Will also accept a single \code{sf} data frame as well.
#' @param move_phase An optional character value indicating a factor variable column in the data that
#' designates different movement phases.
#' @param bm Fit a Brownian Motion model rather than in integrated OU model. Defaults to \code{bm = FALSE}.
#' @param crw_control A named list passed to \code{\link[crawl]{crwMLE}} for optimization.
#' There is one additional parameter \code{lambda}. If added \code{crw_control$lambda}
#' is rate parameter used for the soft equality constraint on class \code{0}, \code{A},
#'  and \code{B} paramters for Argos LS observations in a mixed KF/LS dataset.
#' @param fixPar An alternative to the default set of fixed parameter values. Care should be taken
#' when substituting different values. Make sure you know what you're doing because it can be easily
#' broken
#' @param ... Additional arguments passed to the \code{\link[foreach]{foreach}} function, e.g.,
#' for error handling in the loop.
#' @import dplyr crawl sf foreach progressr
#' @importFrom stats as.formula dexp model.frame model.matrix na.pass qchisq
#' @importFrom utils head
#' @export
#'
cu_crw_argos <- function(data_list, move_phase=NULL, bm=FALSE, crw_control=NULL, fixPar=NULL, ...){
  i <- datetime <- type <- const <- NULL #handle 'no visible binding...'
  progressr::handlers(global = TRUE)
  if(!inherits(data_list,"list")  & inherits(data_list,"sf")){
    data_list <- list(data_list)
  }
  p <- progressr::progressor(length(data_list))
  fits <- foreach(i=1:length(data_list), .packages="sf", ...) %do% {
    dat <- data_list[[i]] %>% dplyr::arrange(datetime)
    all_gps <- all(dat$type%in%c("FastGPS","known"))
    if(all_gps){
      all_ls <- FALSE; all_kf==FALSE; mix_ls_kf <- FALSE
    } else {
      all_ls <- all(dat$type%in%c("Argos_ls","FastGPS","known"))
      all_kf <- all(dat$type%in%c("Argos_kf","FastGPS","known"))
      mix_ls_kf <- !all_ls & !all_kf
    }
    if(!is.null(move_phase)){
      mov.model <- as.formula(paste0("~0+",move_phase))
      mov.mf <-
        model.matrix(mov.model, model.frame(mov.model, dat, na.action = na.pass))
      if (any(is.na(mov.mf)))
        stop("Missing values are not allowed in move_phase!")
      n.mov <- ncol(mov.mf)
    }else{
      n.mov <- 1
      mov.model <- ~1
    }
    theta_mov <- c(rep(8,n.mov),rep(log(-log(0.33)),n.mov))
    fix_mov <- rep(NA, 2*n.mov)
    lc_mov <- c(rep(-Inf,n.mov), rep(log(-log(1-1.0e-4)),n.mov))
    uc_mov <- c(rep(Inf,n.mov), rep(log(-log(1.0e-4)),n.mov))
    if(all_ls){
      err.model <- list(x =  ~0+ln.sd.x+aq0+aqA+aqB)
      fixPar <- c(1,NA,NA,NA,fix_mov)
      constr <- list(
        lower=c(rep(0,3), lc_mov),
        upper=c(rep(Inf,3), uc_mov)
      )
      theta <- c(rep(log(1.2),3),theta_mov)
      prior <- NULL
    } else if(all_kf | all_gps){
      err.model <- list(x =  ~0+ln.sd.x, y = ~0+ln.sd.y, rho= ~error.corr)
      prior <- NULL
      fixPar <- c(1,1,fix_mov)
      constr <- list(
        lower=lc_mov,
        upper=uc_mov
      )
      theta <- theta_mov
      prior <- NULL
    } else if(mix_ls_kf){
      if(is.null(crw_control$lambda)){
        lambda <- 230
      } else {
        lambda <- crw_control$lambda
      }
      err.model <- list(x =  ~0+ln.sd.x+aq0+aqA+aqB, y = ~0+ln.sd.y+aq0+aqA+aqB, rho= ~error.corr)
      prior <- function(par){sum(dexp(abs(par[1:3]-par[4:6]), lambda, log=TRUE))}
      fixPar <- c(1,NA,NA,NA,1,NA,NA,NA,fix_mov)
      constr <- list(
        lower=c(rep(0,6), lc_mov),
        upper=c(rep(Inf,6), uc_mov)
      )
      theta <- c(rep(log(1.2),6),theta_mov)
    } else {
      stop('Unknown location error types!')
    }

    if(bm){
      fixPar <- c(head(fixPar,-n.mov),rep(log(-log(1.0e-4)),n.mov))
      constr$lower <- head(constr$lower,-n.mov)
      constr$upper <- head(constr$upper,-3)
      theta <- head(theta,-3)
    }
    # Fit ctcrw model
    if(is.null(crw_control$initialSANN)){
      initialSANN <- list(maxit=1500, temp=10)
    }else{
      initialSANN <- crw_control$initialSANN
    }
    if(is.null(crw_control$attempts)){
      attempts <- 10
    }else{
      attempts <- crw_control$attempts
    }
    if(is.null(crw_control$control)){
      control <- list(maxit=10000)
    }else{
      control <- crw_control$control
    }
    suppressMessages(
      out <- crawl::crwMLE(
        mov.model = mov.model, err.model = err.model, data = dat, Time.name="datetime",
        fixPar = fixPar, constr = constr, theta = theta,
        control = control, initialSANN = initialSANN,
        prior=prior,
        attempts=attempts, method = "L-BFGS-B")
    )
    p()
    out
  }
  if(length(fits)==1) fits <- fits[[1]]
  return(fits)
}


#' @title Batch CRW Prediction for Multiple Animals
#' @description Uses a list of CRW fitted models and desired prediction times
#' to make location (and velocity) predictions for telemetered animals.
#' @param fit_list A list of CRW fit objects
#' @param predTime A character string describing the desired frequency of prediction,
#' e.g., \code{predTime="1 hour"} or \code{predTime="15 min"}.
#' @param barrier An \code{sf} polygon object representing areas where the animal cannot access.
#' @param vis_graph A visibility graph constructed with the R package \code{pathroutr}, which is used
#' to reroute paths around barriers.
#' @param as_sf Logical. Return an \code{sf} points data frame (\code{TRUE}) or standard \code{crawl} prediction.
#' @param ... Additional arguments passed to the \code{\link[foreach]{foreach}} function, e.g.,
#' for error handling in the loop.
#' @details The R package \code{pathroutr} is necessary for use of the \code{barrier} rerouting.
#' it can be installed with the command
#' \code{install.packages('pathroutr', repos='https://jmlondon.r-universe.dev')}.
#' See 'https://github.com/jmlondon/pathroutr' for a description of use and constructing the
#' viability \code{vis_graph}.
#' @author Devin S. Johnson
#' @export
#' @import sf dplyr crawl foreach progressr
#'
cu_crw_predict <- function(fit_list, predTime=NULL, barrier=NULL, vis_graph=NULL, as_sf=TRUE,...){
  i <- locType <- NULL #handle 'no visible binding...'
  progressr::handlers(global = TRUE)
  route <- !is.null(barrier) & !is.null(vis_graph)
  p <- progressr::progressor(length(fit_list))
  plist <- foreach(i=1:length(fit_list), .packages=c("sf","dplyr"), ...) %do% {
    pred <- crawl::crwPredict(fit_list[[i]], predTime=predTime, return.type="flat")
    if(route){
      if (!requireNamespace("pathroutr", quietly = TRUE)) stop("Please install {pathroutr}: install.packages('pathroutr',repos='https://jmlondon.r-universe.dev')")
      pred <- pred %>% crawl::crw_as_sf(ftype="POINT", locType="p") %>% filter(locType=="p")
      pred <- pred %>% pathroutr::prt_trim(barrier)
      fix <- pathroutr::prt_reroute(pred, barrier, vis_graph, blend=FALSE)
      pred <- pathroutr::prt_update_points(fix, pred)
    }
    if(as_sf & !route) pred <- crw_as_sf(pred,ftype="POINT")
    p()
    pred
  }
  return(plist)
}


#' @title Batch CRW Posterior Path Simulation For Multiple Animals
#' @description Uses a list of CRW fitted models and desired simulation times
#' to make draws from the location (and velocity) posterior distribution for telemetered animals.
#' @param size The number of posterior draws. Defaults to 8 (See Details).
#' @param fit_list A list of CRW fit objects
#' @param predTime A character string describing the desired frequency of prediction,
#' e.g., \code{predTime="1 hour"} or \code{predTime="15 min"}.
#' @param barrier An \code{sf} polygon object representing areas where the animal cannot access.
#' @param vis_graph A visibility graph constructed with the R package \code{pathroutr}, which is used
#' to reroute paths around barriers.
#' @param as_sf Logical. Return an \code{sf} points data frame list (\code{TRUE}) or standard \code{crawl} prediction list
#' @param ... Additional arguments passed to the \code{\link[foreach]{foreach}} function, e.g.,
#' for error handling in the loop.
#' @details The R package \code{pathroutr} is necessary for use of the \code{barrier} rerouting.
#' it can be installed with the command
#' \code{install.packages('pathroutr', repos='https://jmlondon.r-universe.dev')}.
#' See 'https://github.com/jmlondon/pathroutr' for a description of use and constructing the
#' viability \code{vis_graph}.
#' @author Devin S. Johnson
#' @export
#' @import sf dplyr foreach crawl progressr
#'
cu_crw_sample <- function(size=8, fit_list, predTime=NULL, barrier=NULL, vis_graph=NULL, as_sf=TRUE,...){
  i <- j <- NULL #handle 'no visible binding...'
  progressr::handlers(global = TRUE)
  route <- !is.null(barrier) & !is.null(vis_graph)
  p <- progressr::progressor(length(fit_list))
  slist <- foreach(i=1:length(fit_list), .packages=c("sf","dplyr"),...)%do%{
    simObj <- crawl::crwSimulator(fit_list[[i]], parIS = 0, predTime=predTime)
    out <- foreach(j=1:size)%do%{
      samp <- crawl::crwPostIS(simObj, fullPost = FALSE)
      if(route){
        if (! requireNamespace("pathroutr", quietly = TRUE)) stop("Please install {pathroutr}: install.packages('pathroutr',repos='https://jmlondon.r-universe.dev')")
        samp <- samp %>% crawl::crw_as_sf(ftype="POINT", locType="p")
        samp <- samp %>% pathroutr::prt_trim(barrier)
        fix <- pathroutr::prt_reroute(samp, barrier, vis_graph, blend=FALSE)
        samp <- pathroutr::prt_update_points(fix, samp) %>% dplyr::mutate(rep=j)
      }
      if(as_sf & !route) samp <- crw_as_sf(samp,ftype="POINT") %>% dplyr::mutate(rep=j)
      samp
    }
    # if(as_sf) out <- do.call(rbind, out) %>% st_as_sf()
    p()
    out
  }
  return(slist)
}

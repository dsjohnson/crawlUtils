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
  error.corr <- quality <- NULL #handle 'no visible binding...'
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
#' @param bm Fit a Brownian Motion model rather than in integrated OU model. Defaults to \code{bm = FALSE}.
#' @param fixPar An alternative to the default set of fixed parameter values. Care should be taken
#' when substituting different values. Make sure you know what you're doing because it can be easily
#' broken
#' @param ... Additional arguments passed to the \code{\link[foreach]{foreach}} function, e.g.,
#' for error handling in the loop.
#' @import dplyr crawl sf foreach
#' @export
#'
cu_crw_argos <- function(data_list, bm=FALSE, fixPar=NULL, ...){
  i <- datetime <- type <- const <- NULL #handle 'no visible binding...'
  if(!inherits(data_list,"list")  & inherits(data_list,"sf")){
    data_list <- list(data_list)
  }
  fits <- foreach(i=1:length(data_list), .packages="sf", ...) %do% {
    dat <- data_list[[i]] %>% dplyr::arrange(datetime)
    alsg <- all(dat$type%in%c("Argos_ls","FastGPS","known"))
    akfg <- all(dat$type%in%c("Argos_kf","FastGPS","known"))
    if(!alsg & !akfg) {
      pct_alsg <- (dat %>% filter(type=="Argos_ls") %>% nrow()) / nrow(dat)
      pct_akfg <- (dat %>% filter(type=="Argos_kf") %>% nrow()) / nrow(dat)
      argos_pick <- if_else(pct_alsg > pct_akfg, "LS", "KF")
      if (argos_pick == "LS") {
        dat <- dat %>% filter(type %in% c("Argos_ls","FastGPS","known"))
      }
      if (argos_pick == "KF") {
        dat <- dat %>% filter(type %in% c("Argos_kf","FastGPS","known"))
      }
      warning("Animal ", i, " has both LS and KF Argos location types or other unknown types!\n",
              "Keeping ", argos_pick, " becuase ", argos_pick, " represents ",
              "the larger percentage of the observations.")
    }
    if(alsg & akfg) alsg <- FALSE
    if(alsg){
      err.model <- list(x =  ~0+ln.sd.x+aq0+aqA+aqB)
      fixPar <- c(1,NA,NA,NA,NA,NA)
      constr <- list(
        lower=c(rep(0,3), -Inf, log(-log(1-1.0e-4))),
        upper=c(rep(Inf,3), Inf, log(-log(1.0e-4)))
      )
      theta <- c(rep(log(1.2),3),9,log(-log(0.8)))
    } else{
      err.model <- list(x =  ~0+ln.sd.x, y = ~0+ln.sd.y, rho= ~error.corr)
      fixPar <- c(1,1,NA,NA)
      constr <- list(
        lower=c(-Inf, log(-log(1-1.0e-4))),
        upper=c(Inf, log(-log(1.0e-4)))
      )
      theta <- c(9,log(-log(0.8)))
    }
    if(bm){
      fixPar[length(fixPar)] <- log(-log(1.0e-4))
      constr$lower <- constr$lower[-length(constr$lower)]
      constr$upper <- constr$upper[-length(constr$upper)]
      theta <- theta[-length(theta)]
    }
    # Fit ctcrw model
    suppressMessages(
      out <- crawl::crwMLE(
        mov.model = ~1, err.model = err.model, data = dat[897:nrow(dat), ], Time.name="datetime",
        fixPar = fixPar, constr = constr, theta = theta,
        control = list(maxit=10000), initialSANN = list(maxit=1500, temp=10),
        attempts=10, method = "L-BFGS-B")
    )
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
#' @import sf dplyr crawl foreach
#'
cu_crw_predict <- function(fit_list, predTime=NULL, barrier=NULL, vis_graph=NULL, as_sf=TRUE,...){
  i <- locType <- NULL #handle 'no visible binding...'
  route <- !is.null(barrier) & !is.null(vis_graph)
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
#' @import sf dplyr foreach crawl
#'
cu_crw_sample <- function(size=8, fit_list, predTime=NULL, barrier=NULL, vis_graph=NULL, as_sf=TRUE,...){
  i <- j <- NULL #handle 'no visible binding...'
  route <- !is.null(barrier) & !is.null(vis_graph)
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
    out
  }
  return(slist)
}

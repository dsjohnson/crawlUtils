#' @title Predicate function for st_filter
#' @description Predicate function to use with \code{st_filter} such that
#' such that elemets of one spatial object are selected if
#' they are not contained at all in the other. See \link{pkg:sf}{st_within}
#' @param x object of class sf, sfc or sfg
#' @param y object of class sf, sfc or sfg; if missing, x is used
#' @param sparse ogical; should a sparse index list be returned (TRUE) or a dense logical matrix? See \link{pkg:sf}{st_within}.
#' @param prepared ogical; prepare geometry for x, before looping over y? See \link{pkg:sf}{st_within}.
#' @param ... passed on to s2_options
#' @import sf
#' @export
#'
st_not_within <- function(x,y,sparse=TRUE,prepared=TRUE,...){
  !sf::st_within(x,y,sparse,prepared,...)
}


## @title Convert sf output from crawl to ctmm object
## @description Takes \code{crawl::crwPredict} or \code{crawl::crwPostIS}
## objects converted to \code{sf} objects by \code{crawl::crw_as_sf} to data
## objects that can be used within the \code{ctmm} package.
## @param data \code{sf} data object
## @param time.name Character. The name of the column of location times
## @param silent Locical. Should messages from \code{ctmm} be printed. Defaults to
## \code{silent = TRUE}.
## @author Devin S. Johnson
## @import ctmm sf dplyr magrittr
## @export
# crw2ctmm <- function(data, time_name, silent=TRUE){
#   tmp <- data  %>% sf::st_transform(4326) %>% sf::st_coordinates() %>%
#     cbind(data,.) %>% dplyr::select(.data[[time_name]], X, Y, rep) %>% sf::st_drop_geometry() %>%
#     mutate(timestamp=as.character(.data[[time_name]])) %>%
#     rename(
#       `location-long` = X,
#       `location-lat` = Y,
#       `tag-local-identifier` = rep
#     )
#   if(silent){
#     out <- suppressMessages(
#       ctmm::as.telemetry(
#         tmp,
#         projection=sf::st_crs(data)$proj4string)
#     )
#   } else{
#     out <- ctmm::as.telemetry(
#       tmp,
#       projection=sf::st_crs(data)$proj4string
#     )
#   }
#   return(out)
# }


#' @title Label telemetry bout segments
#' @description Label segments of telemetry data such that gaps between observed
#' locations are less than or equal to the interval specified.
#' @param x column of POSIX times for which locations are observed.
#' @param gap specified maximum time interval for which a new bout is started. Defaults
#' to \code{gap=7} days.
#' @param time_unit Unit of time of the gap specification. Defaults to \code{"day"}.
#' @import lubridate units
#' @author Devin S. Johnson
#' @export
#'
cu_get_bouts <- function(x, gap=7, time_unit="day"){
  dt <- diff(x) %>% `units<-`(time_unit)
  time_diff <- c(0, dt)
  bout_id <- (time_diff >= gap) %>% {cumsum(.)+1}
  return(bout_id)
}


#' @title Join crawl prediction or simulation output with a table based on
#' a time interval
#' @description Takes a data set with a POSIX time column named 'datetime' and
#' another data set with 'start' and 'end' columns representing time intervals and
#' merges the two depending whether or not the 'datetime' column is within the interval  of the second.
#' @param x A data frame with a column labeled 'datetime'
#' @param int_tbl A data frame with 'start' and 'end' columns that form non-overlapping intervals as well as at least
#' one other column with interval level data.
#' @author Devin S. Johnson
#' @importFrom fuzzyjoin fuzzy_left_join
#' @export
#'
cu_join_interval_tbl <- function(x, int_tbl){
  x <- fuzzyjoin::fuzzy_left_join(x,seg_tbl,
                                  by=c(datetime="start",datetime="end"),
                                  match_fun = list(`>=`, `<=`)
  )
}


#' @title Location rate statistics
#' @description Calculate location rate statistics such as mean location rate per
#' day or maximum number of locations per day for a telemetry
#' @param x data set containing time of lcoations.
#' @param time_name Character name of the POSIX time column when locations are
#' observed.
#' @param time_unit Time unit of the location summary. Defaults to \code{"day"}.
#' @param stat Function used to summarize location times. Defaults to \code{mean}.
#' @param ... Additional arguments passed to \code{stat} function.
#' @author Devin S. Johnson
#' @import lubridate
#' @export
#'
cu_location_rate <- function(x, time_name, time_unit="day", stat=mean, ...){
  out <- table(lubridate::round_date(x[,time_name], time_unit))
  return(stat(out,...))
}


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
#'
cu_add_argos_cols <- function(x){
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
    x$ln.sd.x <- NA; x$ln.sd.y <- NA; x$error.corr <- 0
  }
  kf_ind <- !(is.na(x$ln.sd.x) | is.na(x$ln.sd.y) | is.na(x$error.corr))
  x <- x %>% dplyr::mutate(
    type = case_when(
      type==type%in%c("Argos","argos","KF") & kf_ind ~ "Argos_kf",
      type==type%in%c("Argos","argos","LS") & !kf_ind ~ "Argos_ls",
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
    error.corr = ifelse(is.na(error.corr), 0, error.corr),
    gq5 = ifelse(quality %in% c("4","5"), 1, 0),
    gq4 = ifelse(quality=="4", 1, 0),
    aq0 = ifelse(quality %in% c("0","A","B"), 1, 0),
    aqA = ifelse(quality %in% c("A","B"), 1, 0),
    aqB = ifelse(quality=="B", 1, 0)
  )
  return(x)
}



#' @title Batch Fitting CTCRW Models for Argos (and FastGPS) Data
#' @description A basic CTCRW model is fitted to a list of data sets where each
#' element in the list represents the data for a single individual or deployment.
#' @param data_list A list of data sets
#' @param bm Fit a Brownian Motion model rather than in integrated OU model. Defaults to \code{bm = FALSE}.
#' @param fixPar An alternative to the default set of fixed parameter values. Care should be taken
#' when substituting different values. Make sure you know what you're doing because it can be easily
#' broken
#' @import dplyr crawl sf progressr foreach
#'
cu_crw_argos <- function(data_list, bm=FALSE, fixPar=NULL){
  p <- progressor(length(data_list))
  fits <- foreach(i=1:length(data_list), .packages="sf") %dorng% {
    dat <- data_list[[i]] %>% dplyr::arrange(datetime)
    alsg <- all(dat$type%in%c("Argos_ls","FastGPS","known"))
    akfg <- all(dat$type%in%c("Argos_kf","FastGPS","known"))
    if(!alsg & !akfg) stop("Animal ", i, " has both LS and KF Argos location types or other unknown types!")
    if(alsg & akfg) alsg <- FALSE
    if(alsg){
      err.model <- list(x =  ~0+ln.sd.x+aq0+aqA+aqB)
      fixPar <- c(1,NA,NA,NA,NA,NA)
      constr <- list(
        lower=c(rep(0,3), -Inf, log(-log(1-1.0e-4))),
        upper=c(rep(Inf,3), Inf, log(-log(1.0e-4)))
      )
      theta <- c(rep(log(1.2),3),9,0.5)
    } else{
      err.model <- list(x =  ~0+ln.sd.x, y = ~0+ln.sd.y, rho= ~error.corr)
      fixPar <- c(1,1,NA,NA)
      constr <- list(
        lower=c(-Inf, log(-log(1-1.0e-4))),
        upper=c(Inf, log(-log(1.0e-4)))
      )
      theta <- c(9,.5)
    }
    if(bm){
      fixPar[length(fixPar)] <- log(-log(1.0e-4))
      constr$lower <- const$lower[-length(const$lower)]
      constr$lower <- const$lower[-length(const$lower)]
      theta <- theta[-length(theta)]
    }
    # Fit ctcrw model
    suppressMessages(
      out <- crawl::crwMLE(
        mov.model = ~1, err.model = err.model, data = dat, Time.name="datetime",
        fixPar = fixPar, constr = constr, theta = theta,
        control = list(maxit=2000), initialSANN = list(maxit=1500, temp=10),
        attempts=10, method = "L-BFGS-B")
    )
    p()
    out
  }
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
#' @details The R package \code{pathroutr} is necessary for use of the \code{barrier} rerouting.
#' it can be installed with the command
#' \code{install.packages('pathroutr', repos='https://jmlondon.r-universe.dev')}.
#' See 'https://github.com/jmlondon/pathroutr' for a description of use and constructing the
#' viability \code{vis_graph}.
#' @author Devin S. Johnson
#' @export
#' @import sf dplyr progressr foreach crawl
#'
cu_batch_predict <- function(fit_list, predTime, barrier=NULL, vis_graph=NULL){
  p <- progressor(length(fit_list))
  plist <- foreach(i=1:length(fit_list), .packages=c("sf","dplyr")) %dorng% {
    pred <- crawl::crwPredict(fit_list[[i]], predTime=predTime, return.type="flat")
    if(!is.null(barrier) & !is.null(vis_graph)){
      if (!requireNamespace("pathroutr", quietly = TRUE)) stop("Please install pathroutr: install.packages('pathroutr',repos='https://jmlondon.r-universe.dev')")
      pred <- pred %>% crawl::crw_as_sf(ftype="POINT", locType="p")
      pred <- pred %>% pathroutr::prt_trim(barrier)
      fix <- pathroutr::prt_reroute(pred, barrier, vis_graph, blend=FALSE)
      pred <- pathroutr::prt_update_points(fix, pred)
    }
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
#' @details The R package \code{pathroutr} is necessary for use of the \code{barrier} rerouting.
#' it can be installed with the command
#' \code{install.packages('pathroutr', repos='https://jmlondon.r-universe.dev')}.
#' See 'https://github.com/jmlondon/pathroutr' for a description of use and constructing the
#' viability \code{vis_graph}.
#' @author Devin S. Johnson
#' @export
#' @import sf dplyr progressr foreach crawl
#'
cu_crw_sample <- function(size=8, fit_list, predTime, barrier=NULL, vis_graph=NULL){
  p <- progressor(length(fit_list))
  slist <- foreach(i=1:length(fit_list), .packages=c("sf","dplyr"))%dorng%{
    simObj <- crawl::crwSimulator(fit_list[[i]], parIS = 0, predTime=predTime)
    out <- foreach(j=1:size)%do%{
      samp <- crawl::crwPostIS(simObj, fullPost = FALSE)
      if(!is.null(barrier) & !is.null(vis_graph)){
        if (! requireNamespace("pathroutr", quietly = TRUE)) stop("Please install pathroutr: install.packages('pathroutr',repos='https://jmlondon.r-universe.dev')")
        samp <- samp %>% crawl::crw_as_sf(ftype="POINT", locType="p")
        samp <- samp %>% pathroutr::prt_trim(barrier)
        fix <- pathroutr::prt_reroute(samp, barrier, vis_graph, blend=FALSE)
        samp <- pathroutr::prt_update_points(fix, samp) %>% dplyr::mutate(rep=j)
      }
      samp
    }
    p()
    out
  }
  return(slist)
}


#' @title Calculate Gaussian Kernel Utilization Distribution
#' @description Calculates a kernel density estimate using telemetry locations
#' over a grid accounting for barriers by normalizing each kernel before addition
#' to the overall sum. The final estimate will always be normalized to sum to one
#' over the estimate grid.
#' @param x an sf or sfc object. The KDE grid will be estimated from the bounding box.
#' @param barrier An sf polygon object that defines an area where use is excluded.
#' @param norm Logical. Should each individual kernel be normalized to account for barriers
#' before addition to the total KDE. Defaults to \code{TRUE} if a barrier is specified.
#' @param ... Any arguments passed to \code{sf::st_make_grid} to create the KDE prediction grid. See \link{pkg:sf}{st_make_grid} for
#' description of additional arguments to make the KDE grid.
#' @import sf
#' @useDynLib crawlUtils, .registration = TRUE
#' @author Devin S. Johnson
#' @export
#'
cu_kde_ud <- function(x, barrier=NULL, norm, ...){
  grid <- sf::st_make_grid(x,...)
}

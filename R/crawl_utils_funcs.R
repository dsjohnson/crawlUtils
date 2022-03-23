#' @title Predicate function for st_filter
#' @description Predicate function to use with \code{st_filter} such that
#' such that elemets of one spatial object are selected if
#' they are not contained at all in the other. See \link{sf::st_within}
#' @param x object of class sf, sfc or sfg
#' @param y object of class sf, sfc or sfg; if missing, x is used
#' @param sparse ogical; should a sparse index list be returned (TRUE) or a dense logical matrix? See \link{sf::st_within}.
#' @param prepared ogical; prepare geometry for x, before looping over y? See \link{sf::st_within}.
#' @param ... passed on to s2_options
#' @import sf
#' @export
not_within <- function(x,y,sparse=TRUE,prepared=TRUE,...){
  !sf::st_within(x,y,sparse,prepared,...)
}


#' @title Convert sf output from crawl to ctmm object
#' @description Takes \code{crawl::crwPredict} or \code{crawl::crwPostIS}
#' objects converted to \code{sf} objects by \code{crawl::crw_as_sf} to data
#' objects that can be used within the \code{ctmm} package.
#' @param data \code{sf} data object
#' @param time.name Character. The name of the column of location times
#' @param silent Locical. Should messages from \code{ctmm} be printed. Defaults to
#' \code{silent = TRUE}.
#' @author Devin S. Johnson
#' @import ctmm sf dplyr magrittr
#' @export
crw2ctmm <- function(data, time_name, silent=TRUE){
  tmp <- data  %>% sf::st_transform(4326) %>% sf::st_coordinates() %>%
    cbind(data,.) %>% dplyr::select(.data[[time_name]], X, Y, rep) %>% sf::st_drop_geometry() %>%
    mutate(timestamp=as.character(.data[[time_name]])) %>%
    rename(
      `location-long` = X,
      `location-lat` = Y,
      `tag-local-identifier` = rep
    )
  if(silent){
    out <- suppressMessages(
      ctmm::as.telemetry(
        tmp,
        projection=sf::st_crs(data)$proj4string)
    )
  } else{
    out <- ctmm::as.telemetry(
      tmp,
      projection=sf::st_crs(data)$proj4string
    )
  }
  return(out)
}


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
get_bouts <- function(x, gap=7, time_unit="day"){
  dt <- diff(x) %>% `units<-`(time_unit)
  time_diff <- c(0, dt)
  bout_id <- (time_diff >= gap) %>% {cumsum(.)+1}
  return(bout_id)
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
location_rate <- function(x, time_name, time_unit="day", stat=mean, ...){
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
#' \code{quality} must be '3','2','1','0','A',or 'B' for \code{type=='Argos'} locations.
#' For other types it is ignored.
#' @author Devin S. Johnson
#' @import dplyr crawl
add_argos_cols <- function(x){
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
  x <- x %>% dplyr::mutate(
    ln.sd.x = dplyr::case_when(
      type == "known" ~ log(10),
      type %in% c("FastGPS","GPS","G","gps") ~ log(50),
      is.na(ln.sd.x) & quality=="3" ~ log(250),
      is.na(ln.sd.x) & quality=="2" ~ log(500),
      is.na(ln.sd.x) & quality %in% c("1","0","A","B") ~ log(1500),
      TRUE ~ ln.sd.x
    ),
    ln.sd.y = case_when(
      type == "known" ~ log(10),
      type %in% c("FastGPS","GPS","G","gps") ~ log(50),
      is.na(ln.sd.y) & quality=="3" ~ log(250),
      is.na(ln.sd.y) & quality=="2" ~ log(500),
      is.na(ln.sd.y) & quality %in% c("1","0","A","B") ~ log(1500),
      TRUE ~ ln.sd.y
    ),
    error.corr = ifelse(is.na(error.corr), 0, error.corr),
    aq0 = ifelse(quality %in% c("0","A","B"), 1, 0),
    aqA = ifelse(quality %in% c("A","B"), 1, 0),
    aqB = ifelse(quality=="B", 1, 0)
  )


}

crwMLE_batch <- function(tdata){
  p <- progressor(nrow(tdata))
  fits <- foreach(i=1:nrow(tdata), .packages="sf") %dorng%{
    temp_dat <- tdata$data[[i]] %>% dplyr::arrange(timestamp)

    # Establish fixed parameter values
    if(temp_dat$loc_type[1] == "Least squares"){
      err.model <- list(
        x =  ~ 0 + ln.sd.x + lc0 + lcA + lcB#,
        # y =  ~ 0 + ln.sd.y + lc0 + lcA + lcB,
        # rho =  ~ error.corr
      )
      fixPar <- c(1,NA,NA,NA,NA,NA)
      constr <- list(
        lower=c(rep(0,3), -Inf, -Inf),
        upper=c(rep(Inf,3), Inf, Inf)
      )
      theta <- c(rep(log(1.2),3),9,0)
    } else {
      err.model <- list(
        x =  ~ 0 + ln.sd.x,
        y =  ~ 0 + ln.sd.y,
        rho =  ~ error.corr)
      fixPar <- c(1,1, NA,NA)
      constr <- list(lower=-Inf, upper=Inf)
      theta <- c(9,0)
    }

    # Fit ctcrw model
    suppressMessages(
      out <- crawl::crwMLE(
        mov.model = ~1, err.model = err.model, data = temp_dat, Time.name="timestamp",
        fixPar = fixPar, constr = constr, theta = theta,
        control = list(maxit=2000), initialSANN = list(maxit=1500, temp=10),
        attempts=10, method = "L-BFGS-B")
    )
    p()
    out
  }
  return(fits)
}



crwPredict_batch <- function(fit_list, predTime, barrier, vis_graph){
  p <- progressor(length(fit_list))
  pred <- foreach(i=1:length(fit_list), .packages=c("sf","dplyr"))%dorng%{
    pred <- crawl::crwPredict(fit_list[[i]], predTime=predTime, return.type="flat") %>%
      crawl::crw_as_sf(ftype="POINT", locType="p") %>%
      # mutate(
      #   tmp_seg1 = zoo::na.locf(seg_id),
      #   tmp_seg2 = zoo::na.locf(seg_id, fromLast=TRUE),
      #   seg_id = ifelse(tmp_seg1!=tmp_seg2, NA, tmp_seg1)
      # ) %>% select(-tmp_seg1, -tmp_seg2) %>%
      pathroutr::prt_trim(iso20)
    fix <- pathroutr::prt_reroute(pred, barrier, vis_graph, blend=FALSE)
    pred <- pathroutr::prt_update_points(fix, pred)
    p()
    pred
  }
  return(pred)
}


crwSample_batch <- function(size, fit_list, predTime, barrier, vis_graph){
  p <- progressor(length(fit_list))
  slist <- foreach(i=1:length(fit_list), .packages=c("sf","dplyr"))%dorng%{
    simObj <- crawl::crwSimulator(fit_list[[i]], parIS = 0, predTime="15 min")
    out <- foreach(j=1:size)%do%{
      samp <- crawl::crwPostIS(simObj, fullPost = FALSE) %>%
        crawl::crw_as_sf(ftype="POINT", locType="p") %>%
        pathroutr::prt_trim(iso20)
      samp_fix <- pathroutr::prt_reroute(samp, iso20, vis_graph, blend=FALSE)
      samp <- pathroutr::prt_update_points(samp_fix, samp) %>% dplyr::mutate(rep=j)
      samp
    }
    p()
    out
  }
  return(slist)
}

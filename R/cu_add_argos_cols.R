#' @title Add columns for modeling FastGPS and ARGOS error structure in the same
#' telemetry deployment
#' @description Columns are added to the telemetry data set so that multiple data types
#' can be used simultaneously within the same animal: FastGPS, Argos least-squares, and
#' Argos Kalman filter.
#' @param x Data frame containing location telemetry data and Argos quality
#' information. See 'Details' for a description of the necessary data column names.
#' @param units Units for movement and location error models. One of `"meter"`
#' (`"metre"`) or `"kilometer"` (`"kilometre"`). If the data are projected, it will automatically use the units
#' of the locations and ignore this argument.
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
cu_add_argos_cols <- function(x, units="meter"){
  if (inherits(x,"sf")) {
    proj_unit <- sf::st_crs(x, parameters = TRUE)$units_gdal
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

  error.corr <- quality <- ln.sd.x <- ln.sd.y <- error_area <- NULL #handle 'no visible binding...'
  col_nms <- colnames(x)
  if('error_semi_major_axis' %in% col_nms &
     'error_semi_minor_axis' %in% col_nms &
     'error_ellipse_orientation' %in% col_nms){
    x <- cbind(x,
               cu_argosDiag2Cov(
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
        type == "known" ~ log(units*10),
        type=="FastGPS" & quality=="4" ~ log(1163*units/2),
        type=="FastGPS" & quality=="5" ~ log(169*units/2),
        type=="FastGPS" & quality=="6" ~ log(71*units/2),
        type=="FastGPS" & quality=="7" ~ log(43*units/2),
        type=="FastGPS" & quality=="8" ~ log(34*units/2),
        type=="FastGPS" & quality=="9" ~ log(28*units/2),
        type=="FastGPS" & quality=="10" ~ log(24*units/2),
        type=="FastGPS" & quality=="11" ~ log(19*units/2),
        type=="Argos_ls" & quality=="3" ~ log(250*units),
        type=="Argos_ls" & quality=="2" ~ log(500*units),
        type=="Argos_ls" & quality=="1" ~ log(1500*units),
        type=="Argos_ls" & quality=="0" ~ log(2.25*1500*units),
        type=="Argos_ls" & quality=="A" ~ log(3.98*1500*units),
        type=="Argos_ls" & quality=="B" ~ log(7.37*1500*units),
        TRUE ~ ln.sd.x
      ),
      ln.sd.y = dplyr::case_when(
        type == "known" ~ log(units*10),
        type=="FastGPS" & quality=="4" ~ log(1163*units/2),
        type=="FastGPS" & quality=="5" ~ log(169*units/2),
        type=="FastGPS" & quality=="6" ~ log(71*units/2),
        type=="FastGPS" & quality=="7" ~ log(43*units/2),
        type=="FastGPS" & quality=="8" ~ log(34*units/2),
        type=="FastGPS" & quality=="9" ~ log(28*units/2),
        type=="FastGPS" & quality=="10" ~ log(24*units/2),
        type=="FastGPS" & quality=="11" ~ log(19*units/2),
        type=="Argos_ls" & quality=="3" ~ log(250*units),
        type=="Argos_ls" & quality=="2" ~ log(500*units),
        type=="Argos_ls" & quality=="1" ~ log(1500*units),
        type=="Argos_ls" & quality=="0" ~ log(2.5*1500*units),
        type=="Argos_ls" & quality=="A" ~ log(3.67*1500*units),
        type=="Argos_ls" & quality=="B" ~ log(5.42*1500*units),
        TRUE ~ ln.sd.y
      ),
      error.corr = ifelse(is.na(.data$error.corr), 0, .data$error.corr),
      # gq5 = ifelse(quality %in% c("5"), 1, 0),
      # gq4 = ifelse(quality=="4", 1, 0),
      # aq0 = ifelse(quality %in% c("0","A","B"), 1, 0),
      # aqA = ifelse(quality %in% c("A","B"), 1, 0),
      # aqB = ifelse(quality=="B", 1, 0),
      low_qual_argos = ifelse(quality %in% c("0","A","B"), 1, 0),
      low_qual_gps = ifelse(quality=="4", 1, 0),
      error_area = pi*qchisq(0.95, 2)*exp(0.5*(ln.sd.x+ln.sd.y))*sqrt(1-error.corr^2),
      # error_area = ifelse(quality=="0", 1.1*error_area, error_area),
      # error_area = ifelse(quality=="A", 1.2*error_area, error_area),
      # error_area = ifelse(quality=="B", 1.3*error_area, error_area)
    )
  return(x)
}


get_ls_error_terms <- function(data){
  mod <- NULL
  par <- 0
  if(any(data$quality=='0')){
    mod <- " + aq0"
    par <- par+1
  }
  if(any(data$quality=='A')){
    mod <- paste0(mod," + aqA")
    par <- par+1
  }
  if(any(data$quality=='B')){
    mod <- paste0(mod," + aqB")
    par <- par+1
  }
  return(list(mod=mod, par=par))
}


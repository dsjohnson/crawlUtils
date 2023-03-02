#' @title Add cluster information from the {GPSeq_clus} package.
#' @description The function is simply a wrapper for the \code{\link[GPSeqClus]{GPSeq_clus}} function. All of the
#' arguments are the same except for the \code{dat} argument. The \code{dat} argument
#' must be an object output from \code{\link[crawlUtils]{cu_crw_predict}} or \code{\link[crawlUtils]{cu_crw_sample}} with
#' \code{as_sf=TRUE}. In addition, here \code{show_plots=c(FALSE, "mean")} as the default. As the function executes
#' there are a number of progress bars, they can be ignored.
#' @param dat Output from \code{cu_crw_predict} or \code{cu_crw_sample}.
#' @param search_radius_m Search radius (meters) from cluster centroid when building clusters.
#' @param window_days Temporal window (days) to search for new locations from the most recent location in a cluster
#' @param clus_min_locs Minimum number of locations required to form a cluster. Default is 2.
#' @param centroid_calc Method for recalculating centroids when actively building clusters - e.g., "median" or "mean" (default). Not to be confused with
#'                      plotting the "mean" or "median" centroid once a cluster has been built.
#' @param show_plots Vector of TRUE/FALSE for plotting followed by plotting argument for the "median" or "mean" centroid - e.g., c(TRUE, "mean") (default)
#' @param scale_plot_clus When plotting, scale cluster markers based on number of locations (TRUE/FALSE).
#' @param store_plots When plotting, also assign map outputs to global environment (TRUE/FALSE).
#' @param season_breaks_jul Ascending numeric vector of julian days (0-365) used to classify by season/parturition/hunting seasons etc.
#'                          e.g., c(121, 274, 305) result may be: 1 Nov - 30 Apr (winter = 0), 1 May - 31 Aug (summer = 1), 1 Oct - 31 Oct (hunting season = 2)
#' @param daylight_hrs Manually set start and stop hours (0-24) to classify day and night locations. - e.g. c(6,18) would classify 6AM - 6PM as daylight hrs.
#'                     NA (default) uses 'suncalc' package to convert cluster location and time to be classified based on specific specific sunrise and sunset times.
#' @references Clapp, J. G., Holbrook, J. D., & Thompson, D. J. (2021). GPSeqClus:
#' An R package for sequential clustering of animal location data for model
#' building, model application and field site investigations. Methods in Ecology
#'  and Evolution, 12(5), 787-793.
#' @import sf dplyr
#' @importFrom GPSeqClus GPSeq_clus
#' @importFrom sfheaders sf_to_df
#' @importFrom rlang .data
#' @author Devin S. Johnson
#' @export
cu_GPSeq_clus <- function(dat, search_radius_m, window_days, clus_min_locs = 2,
                          centroid_calc = "mean", show_plots = c(FALSE, "mean"), scale_plot_clus = TRUE,
                          store_plots = FALSE, season_breaks_jul = NA, daylight_hrs = NA){

  data <- AID <- TelemDate <- datetime <- y <- point_id <- sfg_id <- . <- .data <- x <- NULL

  ## Function to convert crw objects to GPSeqClus objects
  crw_to_GPSeq <- function(dat){
    if("x" %in% colnames(dat)) dat <- dplyr::rename(dat, x..1=.data[["x"]])
    if("y" %in% colnames(dat)) dat <- dplyr::rename(dat, y..1=.data[["y"]])
    dat <- st_transform(dat, 4326) |> sfheaders::sf_to_df(fill=TRUE) |> dplyr::select(-sfg_id, -point_id) |>
      dplyr::rename(Long=x, Lat=y, TelemDate=datetime)
    if("x..1" %in% colnames(dat)) dat <- dplyr::rename(dat, x=.data[["x..1"]])
    if("y..1" %in% colnames(dat)) dat <- dplyr::rename(dat, y=.data[["y..1"]])
    dat$AID <- 1
    return(dat)
  }

  ## function to execute clustering
  ex_GPSeq_clus <- function(dat,search_radius_m=search_radius_m, window_days=window_days, clus_min_locs=clus_min_locs,
                            centroid_calc=centroid_calc, show_plots=show_plots, scale_plot_clus=scale_plot_clus,
                            store_plots=store_plots, season_breaks_jul=season_breaks_jul, daylight_hrs=daylight_hrs){
    x_type <- attr(dat, "crw_type")
    x_crs <- st_crs(dat)
    dat <- crw_to_GPSeq(dat)
    clus <- try(suppressWarnings(
      GPSeqClus::GPSeq_clus(dat=dat,search_radius_m=search_radius_m, window_days=window_days, clus_min_locs=clus_min_locs,
                            centroid_calc=centroid_calc, show_plots=show_plots, scale_plot_clus=scale_plot_clus,
                            store_plots=store_plots, season_breaks_jul=season_breaks_jul, daylight_hrs=daylight_hrs)
    ), silent=TRUE)
    if(!inherits(clus,"try-error")){
      clus <- clus[[1]] |> rename(datetime=TelemDate) |> st_as_sf(coords=c("Long","Lat"), crs=4326) |>
        st_transform(x_crs)
    } else {
      clus <- dat
      clus <- clus |> rename(datetime=TelemDate) |> st_as_sf(coords=c("Long","Lat"), crs=4326) |>
        st_transform(x_crs)
      warning("Clustering failed, returning original data. Try different 'search_radius_m' values")
    }
    clus <- select(clus, -AID)
    attr(clus, "crw_type") <- x_type
    return(clus)
  }

  ## Main portion of function
  x_type <- attr(dat, "crw_type")
  if(inherits(dat, "list")){
    if(!all(sapply(dat, attr, "crw_type")=="crwIS_sf")) stop("The 'dat' argument is not the correct form!")
    x_type <- "crwIS_sf_list"
  }
  if(! x_type %in% c("crwIS_sf_list","crwIS_sf","crwPredict_sf")) stop("The 'dat' argument is not the correct form!")
  if(x_type %in% c("crwIS_sf","crwPredict_sf")){
    dat <- ex_GPSeq_clus(dat,  search_radius_m=search_radius_m, window_days=window_days,
                       clus_min_locs=clus_min_locs,centroid_calc=centroid_calc,
                       show_plots=show_plots, scale_plot_clus=scale_plot_clus,
                       store_plots=store_plots, season_breaks_jul=season_breaks_jul,
                       daylight_hrs=daylight_hrs)
  } else{
    for(i in 1:length(dat)){
      dat[[i]] <- ex_GPSeq_clus(dat[[i]],  search_radius_m=search_radius_m, window_days=window_days,
                              clus_min_locs=clus_min_locs,centroid_calc=centroid_calc,
                              show_plots=show_plots, scale_plot_clus=scale_plot_clus,
                              store_plots=store_plots, season_breaks_jul=season_breaks_jul,
                              daylight_hrs=daylight_hrs)
    }
    attr(dat, "crw_type") <- "crwIS_sf_list"
  }
  return(dat)
}

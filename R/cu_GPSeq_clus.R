#' @title Add cluster information from the {GPSeq_clus} package.
#' @description The function is simply a wrapper for the \code{\link[GPSeqClus]{GPSeq_clus}} function. All of the
#' arguments are the same except for the \code{dat} argument (here \code{x}). The \code{x} argument
#' must be an object output from \code{\link[crawlUtils]{cu_crw_predict}} or \code{\link[crawlUtils]{cu_crw_sample}} with
#' \code{as_sf=TRUE}. In addition, here \code{show_plots=c(FALSE, "mean")} as the default. As the function executes
#' there are a number of progress bars, they can be ignored.
#' @param x Output from \code{cu_crw_predict} or \code{cu_crw_sample}.
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
cu_GPSeq_clus <- function(x, search_radius_m, window_days, clus_min_locs = 2,
                          centroid_calc = "mean", show_plots = c(FALSE, "mean"), scale_plot_clus = TRUE,
                          store_plots = FALSE, season_breaks_jul = NA, daylight_hrs = NA){
  data <- AID <- TelemDate <- datetime <- y <- point_id <- sfg_id <- . <- .data <- NULL
  x_type <- attr(x, "crw_type")
  if(inherits(x, "list")){
    if(!all(sapply(x, attr, "crw_type")=="crwIS_sf")) stop("The 'x' argument is not the correct form!")
    x <- x %>%  do.call(rbind,.)
    attr(x, "crw_type") <- "crwIS_sf_list"
    x_type <- "crwIS_sf_list"
  }
  if(! x_type %in% c("crwIS_sf_list","crwIS_sf","crwPredict_sf")) stop("The 'x' argument is not the correct form!")
  x_crs <- st_crs(x)
  if("x" %in% colnames(x)) x <- rename(x, x..1=.data[["x"]])
  if("y" %in% colnames(x)) x <- rename(x, y..1=.data[["y"]])
  x <- st_transform(x, 4326) |> sfheaders::sf_to_df(fill=TRUE) |> select(-sfg_id, -point_id) |>
    rename(Long=x, Lat=y, TelemDate=datetime)
  if("x..1" %in% colnames(x)) x <- rename(x, x=.data[["x..1"]])
  if("y..1" %in% colnames(x)) x <- rename(x, y=.data[["y..1"]])
  if(x_type=="crwIS_sf_list"){
    x <- mutate(x, AID=rep)
  } else {
    x <- mutate(x, AID=1)
  }
  clus <- suppressWarnings(
    GPSeqClus::GPSeq_clus(dat=x,  search_radius_m=search_radius_m, window_days=window_days, clus_min_locs=clus_min_locs,
                           centroid_calc=centroid_calc, show_plots=show_plots, scale_plot_clus=scale_plot_clus,
                           store_plots=store_plots, season_breaks_jul=season_breaks_jul, daylight_hrs=daylight_hrs)
  )
  clus <- clus[[1]] |> rename(datetime=TelemDate) |> st_as_sf(coords=c("Long","Lat"), crs=4326) |>
    st_transform(x_crs)
  split_f <- clus$AID
  clus <- select(clus, -AID)
  if(length(unique(split_f))>1){
    clus <- split(clus, f=split_f)
    for(i in 1:length(clus)) attr(clus[[i]],"crw_type") <- "crwIS_sf"
  }
  attr(clus, "crw_type") <- x_type
  return(clus)
}

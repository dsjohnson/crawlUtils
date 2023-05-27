#' @title Convert Wildlife Computers data imported with
#' `\link[crawlUtils]{cu_read_wc_dirs}` to a `telemetry` object from the
#' `ctmm` package.
#' @param locs_df A data frame output by the function `\link[crawlUtils]{cu_read_wc_dirs}`
#' @param ... Additional arguments to be passed to `\link{ctmm}{as.temeletry}`
#' @author Josh M. London, Devin S. Johnson
#' @import dplyr
#' @importFrom tibble tibble
#' @importFrom ctmm as.telemetry uere<-
#' @importFrom tidyr nest
#' @export
cu_as_telem <- function(locs_df,...) {
  type <- deploy_id <- datetime <- longitude <- latitude <- error_ellipse_orientation <-
    error_semi_minor_axis <- error_semi_major_axis <- quality <- data <- NULL
  # if(is.null(projection)){
  #   locs_sf <- st_as_sf(locs_df, coords=c("longitude","latitude"), crs=4326)
  #   prj <- crsuggest::suggest_crs(locs_sf)$crs_proj4[[1]]
  #   message("Using projection: ", prj)
  #   rm(locs_sf)
  # } else if(projection=="ctmm_default"){
  #   warning("Using {ctmm} default projections. This will result in a potentially different projection for each animal!")
  #   prj <- NULL
  # } else{
  #   prj <- projection
  # }
  # separate fastloc and argos
  locs_f <- locs_df |> filter(type%in%c("FastGPS","known")) |> rm_dup()
  locs_a <- locs_df |> filter(type=="Argos") |> rm_dup()
  rm(locs_df)

  # rename for movebank conventions and convert
  locs_a <- locs_a |>
    rename(
      individual.local.identifier = deploy_id,
      timestamp = datetime,
      location.long = longitude,
      location.lat = latitude,
      Argos.orientation = error_ellipse_orientation,
      Argos.semi.minor = error_semi_minor_axis,
      Argos.semi.major = error_semi_major_axis
    ) %>% mutate(
      Argos.location.class = quality,
      quality = as.character(quality)
    )
  locs_a <- ctmm::as.telemetry(object = locs_a, ...)
  locs_a <- tibble(deploy_id=names(locs_a), telem=locs_a)

  locs_f <- locs_f |>
    rename(
      individual.local.identifier = deploy_id,
      timestamp = datetime,
      location.long = longitude,
      location.lat = latitude
    ) %>% mutate(
      HDOP = dplyr::case_when(
        type == "known" ~ sqrt(2),
        type=="FastGPS" & quality=="4" ~ sqrt(2)*(1163)/20,
        type=="FastGPS" & quality=="5" ~ sqrt(2)*(169)/20,
        type=="FastGPS" & quality=="6" ~ sqrt(2)*(71)/20,
        type=="FastGPS" & quality=="7" ~ sqrt(2)*(43)/20,
        type=="FastGPS" & quality=="8" ~ sqrt(2)*(34)/20,
        type=="FastGPS" & quality=="9" ~ sqrt(2)*(28)/20,
        type=="FastGPS" & quality=="10" ~ sqrt(2)*(24)/20,
        type=="FastGPS" & quality=="11" ~ sqrt(2),
        TRUE ~ Inf
      ),
      quality = as.character(quality)
    )
  locs_f <- ctmm::as.telemetry(object = locs_f, ...)
  uere(locs_f) <- 20
  locs_f <- tibble(deploy_id=names(locs_f), telem=locs_f)


  locs_df <- bind_rows(locs_a, locs_f) |> group_by(deploy_id) |> nest()
  locs_df <- locs_df |> rowwise() |> mutate(
    data = list(data$telem |> ctmm::tbind())
  )

  names(locs_df$data) <- locs_df$deploy_id

  return(locs_df$data)
}

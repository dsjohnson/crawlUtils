#' @title Movement bout detection
#' @description Creates a data table that indicates the times of different
#' bouts of movement. This method uses changes in the overall dispersion rate of the animal
#' from the 'base' time to detect changes in overall movement from small scale local movement
#' to large scale migration.
#' @param data An \code{sf} data set of locations with times noted by the 'datetime' column
#' @param min_disp The minimum dispersion rate to be considered a migration interval,
#' e.g. 10 for a 10km dispersion minimum.
#' @param  migr_speed_cut The minimum speed to categorize a movement bout as "migratory" Defaults to 1 km/h.
#' @param min_bout_len Minimum length for each bout, defaults to 3
#' @param Gmax Maximum number of movement clusters, defaults to 3.
#' @param grid_res The temporal resolution at which migrations are detected. e.g., "day" (default) implies migration
#' start and end is detected on a daily resolution.
#' @param base The location at which dispersion is measured. Can be one of \code{"first"} (first location),
#' \code{"last"} (final location), or some other \code{sf::sfc} point location.
#' @export
#' @author Devin S. Johnson
#' @import units lubridate dplyr sf
#' @importFrom stats coef dist predict vcov coefficients lm
#' @importFrom mclust Mclust
#'
cu_bout_det_mbc <- function(data, min_disp, migr_speed_cut = 1, min_bout_len=3, Gmax=3,
                             grid_res="day", base="first"){
  travel <- bout <- datetime <- speed <- disp_rate <- NULL
  if(base=="first"){
    base <- data[1,] %>% st_geometry()
  } else if(base=="last"){
    base <- data[nrow(data),] %>% st_geometry()
  } else if(!inherits(base,c("sf","sfc"))){
    stop("Unrecognized 'base' argument!")
  }
  fit <- cu_crw_argos(data, bm=TRUE)
  predTime <- c(data$datetime[1],
                seq(ceiling_date(data$datetime[1], grid_res), floor_date(tail(data$datetime,1), grid_res), by=grid_res),
                tail(data$datetime,1)
                )
  pred <- cu_crw_predict(fit, predTime=predTime)
  pred <- pred[pred$locType=="p",]
  pred$speed <-  set_units(pred$speed, "m/h") %>% set_units("km/h")

  ddd <- data.frame(
    datetime = pred$datetime,
    dist = st_distance(pred, pred[1,]) %>% units::set_units("km") %>% as.vector()
  )
  ddd$time <- with(ddd,as.numeric(datetime)/as.numeric(duration(1, grid_res)))

  # cluster attempt ######
  if(all(ddd$dist<=min_disp)){
    mov_type <- rep(1, nrow(ddd))
  } else{
    clust_fit <- Mclust(cbind(pred$speed, pred$TimeNum, ddd$dist), G=1:Gmax)
    mov_type <- clust_fit$classification
  }
  if(length(unique(mov_type))==1){
    avg_speed <- set_units(mean(pred$speed), "km/h")
  } else{
    avg_speed <- as.vector(coefficients(lm(pred$speed ~ 0 + as.factor(mov_type)))) %>%
      set_units("km/h")

  }
  stat_group <- which(avg_speed <= set_units(migr_speed_cut, "km/h"))
  mov_type <- as.numeric(!mov_type %in% stat_group)
  x <- rle(mov_type)
  ###

  x$values[x$lengths<=min_bout_len & x$values==1] <- 0
  x <- rle(inverse.rle(x))
  x$values[x$lengths<=min_bout_len & x$values==0] <- 1
  x <- rle(inverse.rle(x))
  x$values <- cumsum(x$values)*x$values

  ddd$travel <- as.integer(inverse.rle(x))
  ddd$speed <- pred$speed
  ddd$bout <-with(rle(ddd$travel), rep(seq_along(values), lengths))
  summ <- group_by(ddd, travel, bout) %>%
    summarize(
      start = floor_date(min(datetime)),
      # end = as.Date(max(datetime)),
      avg_speed = mean(speed, na.rm=TRUE),
      .groups="drop"
    ) %>% arrange(bout)
  summ <- rbind(summ, data.frame(travel=NA, bout=NA, start=floor_date(max(ddd$datetime)) + duration(grid_res), avg_speed=NA))

  attr(summ, "base") <- base

  return(summ)
}

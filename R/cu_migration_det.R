#' @title Migration detection
#' @description Creates a data table that indicates the times of different
#' bouts of movement. This method uses changes in the overall dispersion rate of the animal
#' from the 'base' time to detect changes in overall movement from small scale local movement
#' to large scale migration.
#' @param data An \code{sf} data set of locations with times noted by the 'datetime' column
#' @param min_disp The minimum dispersion rate to be considered a migration interval,
#' e.g. 10 for a 10km dispersion minimum.
#' @param  migr_disp_cut The minimum per day dispersal rate to categorize a movement bout as "migratory" Defaults to 1.
#' @param min_bout_len The mimimum length of time that a migration or non-migration event
#' will take, e.g., 7 implies a minimum of 7 time intervals for a bout.
#' @param grid_res The temporal resolution at which migrations are detected. e.g., "day" (default) implies migration
#' start and end is detected on a daily resolution.
#' @param base The location at which dispersion is measured. Can be one of \code{"first"} (first location),
#' \code{"last"} (final location), or some other \code{sf::sfc} point location.
#' @param max_k The maximum degrees of freedom used by \code{mgcv::gam} to model dispersion and estimate
#' the derivative of the dispersion function.
#' @export
#' @author Devin S. Johnson
#' @import units lubridate dplyr sf
#' @rawNamespace import(mgcv, except = mvn)
#' @importFrom stats coef dist predict vcov coefficients lm
#' @importFrom mclust Mclust
#'
cu_migration_det <- function(data, min_disp, migr_disp_cut = 1, min_bout_len=3,
                             grid_res="day", base="first", max_k=100){
  travel_evt <- bout <- datetime <- disp_rate <- NULL
  if(base=="first"){
    base <- data[1,] %>% st_geometry()
  } else if(base=="last"){
    base <- data[nrow(data),] %>% st_geometry()
  } else if(!inherits(base,c("sf","sfc"))){
    stop("Unrecognized 'base' argument!")
  }
  if(nrow(data)<= (max_k+1) ){
    k <- floor(nrow(data)/2)
  } else{
    k <- as.integer(max_k)
  }
  ddd <- data.frame(
    datetime = data$datetime,
    dist = st_distance(data, base) %>% units::set_units("km") %>% as.vector()
  )
  ddd$time <- with(ddd,
                   as.numeric(datetime)/as.numeric(duration(1, grid_res))
  )
  suppressWarnings(dfit <- mgcv::gam(dist~s(time, k=k, bs='ad'), data=ddd, method='REML'))
  grid1 <- floor_date(min(ddd$datetime), grid_res)
  grid2 <- ceiling_date(max(ddd$datetime), grid_res)
  newdata <- data.frame(datetime = seq(grid1, grid2, grid_res))
  newdata$time <- with(newdata,
                       as.numeric(datetime)/as.numeric(duration(1, grid_res))
  )
  X <- predict(dfit, newdata=newdata, type="lpmatrix")
  der <- abs(diff(X%*%coef(dfit)))
  # D <- -1*diag(nrow(newdata))
  # for(i in 1:(nrow(D)-1)){
  #   D[i,i+1] <- 1
  # }
  # D <- D[-nrow(D),]
  # V <- sqrt(diag(D%*%X%*%vcov(dfit, unconditional=T)%*%t(X)%*%t(D)))
  # derup <- der + 1.96*V
  # derlo <- der - 1.96*V
  # sig <- (derup*derlo >0)*(abs(der)>=min_disp)
  # x <- rle(as.numeric(sig))

  # cluster attempt ######
  if(all(ddd$dist<=min_disp)){
    mov_type <- rep(1, length(der))
  } else{
    clust_fit <- Mclust(der, G=1:4)
    mov_type <- clust_fit$classification
  }
  if(length(unique(mov_type))==1){
    avg_disp <- mean(der)
  } else{
    avg_disp <- coefficients(lm(der ~ 0 + as.factor(mov_type)))
  }
  stat_group <- which(avg_disp <= migr_disp_cut)
  mov_type <- as.numeric(!mov_type %in% stat_group)
  x <- rle(mov_type)
  ###

  x$values[x$lengths<=min_bout_len & x$values==1] <- 0
  x <- rle(inverse.rle(x))
  x$values[x$lengths<=min_bout_len & x$values==0] <- 1
  x <- rle(inverse.rle(x))
  # top_runs <- sort(x$lengths[x$values==1], decreasing=TRUE)[1:max_num_mig]
  # trl <- length(x$lengths[x$values==1 & x$lengths%in%top_runs])
  # if(trl>max_num_mig) warning("There were multiple migration events with the same length")
  # x$values[x$values==1 & x$length<min(top_runs)] <- 0
  x$values <- cumsum(x$values)*x$values
  newdata$travel_evt <- c(inverse.rle(x),tail(x$values,1))
  newdata$disp_rate <- c(der, NA)
  newdata$bout <-with(rle(newdata$travel_evt), rep(seq_along(values), lengths))
  summ <- group_by(newdata, travel_evt, bout) %>%
    summarize(
      start = as.Date(min(datetime)),
      end = as.Date(max(datetime)),
      avg_disp_rate = mean(disp_rate, na.rm=TRUE),
      .groups="drop"
    ) %>% arrange(bout)
  summ$end[1:(nrow(summ)-1)] <- summ$end[1:(nrow(summ)-1)]+duration(grid_res)
  summ$avg_disp_rate <- units::set_units(summ$avg_disp_rate, paste0("km/",grid_res), mode='standard')
  summ$travel_evt <- factor(summ$travel_evt)
  summ$bout <- factor(summ$bout)

  attr(summ, "base") <- base

  return(summ)
}

#' @title Read individual telemetry data from Wildlife Computers portal directories
#' @description Read and combine data downloaded from Wildlife Computers portal
#' into individual directories.
#' @param x Directory containing the individual telemetry data directories.
#' @export
#' @author Devin S. Johnson
#' @import lubridate dplyr
#' @importFrom janitor clean_names
#' @importFrom readr read_csv

cu_read_wc_dirs <- function(x){
  #
  Latitude <- Longitude <- NULL
  # Determine which file to load for each animal:
  dirs <- list.dirs(x)[-1]
  nms1 <- paste0(list.dirs(x, full.names=FALSE)[-1],
                 "-1-Locations.csv"
  )
  loc_file1 <- paste(dirs, nms1, sep="/")
  nms2 <- paste0(list.dirs(x, full.names=FALSE)[-1],
                 "-Locations.csv"
  )
  loc_file2 <- paste(dirs, nms2, sep="/")

  # Container for file names
  loc_file <- ifelse(file.exists(loc_file1), loc_file1, loc_file2)

  # Read in data and combine into single table
  # There are 2 animals with no location data
  locs <- NULL
  for(i in 1:length(loc_file)){
    if(!file.exists(loc_file[i])) next
    id_data <- readr::read_csv(loc_file[i], show_col_types=FALSE) %>%
      filter(!is.na(Latitude), !is.na(Longitude))
    time <- parse_date_time(id_data$Date,"%H:%M:%S %d-%b-%Y", quiet=TRUE)
    if(all(is.na(time))){
      time <- parse_date_time(id_data$Date,"%m/%d/%y %H:%M", quiet=TRUE)
    }
    if(all(is.na(time))){
      time <- parse_date_time(id_data$Date,"%m/%d/%y %H:%M:%S", quiet=TRUE)
    }
    if(all(is.na(time))){
      time <- parse_date_time(id_data$Date,"%m/%d/%Y %H:%M:%S", quiet=TRUE)
    }
    if(all(is.na(time))){
      stop("Date format is unrecognizable.")
    }
    if(any(is.na(time))){
      warning("Some dates were not converted to POSIX format. Look for NAs in datetime.")
    }
    id_data$datetime <- time
    id_data <- select(id_data, -Date)
    id_data <- janitor::clean_names(id_data)
    locs <- bind_rows(locs, id_data)
  }
return(locs)
}
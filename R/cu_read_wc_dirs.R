#' @title Read individual telemetry data from Wildlife Computers portal directories
#' @description Read and combine data downloaded from Wildlife Computers portal
#' into individual directories.
#' @param x Directory containing the individual telemetry data directories.
#' @param remove_duplicates Locgical. Should observations with duplicated times be removed? The observation
#' with the highest quality will be retained.
#' @export
#' @author Devin S. Johnson, Josh M. London
#' @import lubridate dplyr
#' @importFrom janitor clean_names
#' @importFrom readr read_csv

cu_read_wc_dirs <- function(x, remove_duplicates=TRUE){
  #
  Latitude <- Longitude <- quality <- NULL
  # Determine which file to load for each animal:
  dirs <- list.dirs(x)[-1]

  # Read in data and combine into single table
  # There are 2 animals with no location data
  locs <- NULL
  for(i in 1:length(dirs)){
    loc_paths <- list.files(dirs[[i]], pattern="*-Locations.csv")
    if(length(loc_paths)==1){
      loc_file <- loc_paths[[1]]
    } else{
      loc_paths <- strsplit(loc_paths,"-")
      loc_paths <- loc_paths[sapply(loc_paths, "length")>2]
      run <- as.numeric(sapply(loc_paths, \(x) x[[2]]))
      loc_file <- paste(loc_paths[run==max(run)][[1]], collapse="-")
    }
    loc_file <- paste0(dirs[[i]],"/",loc_file)
    id_data <- readr::read_csv(loc_file, show_col_types=FALSE) %>%
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
  locs <- locs |> mutate(
    quality = factor(quality, levels=c(as.character(11:0),"A","B","Z"))
  )

  if(remove_duplicates) locs <- rm_dup(locs)

return(locs)
}




#' @title Migration detection
#' @description Uses daily displacement from the initial location to
#' determine stat and stop times of migration intervals.
#' @param data
#' @param expected
#' @param min
#' @param from
#' @param k
#' @export
#' @importFrom mgcv gam
#' @importFrom lubridate floor_date ceiling_date
#'

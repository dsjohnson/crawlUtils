#' @title Calculate Gaussian Kernel Utilization Distribution
#' @description Calculates a kernel density estimate using telemetry locations
#' over a grid accounting for barriers by normalizing each kernel before addition
#' to the overall sum. The final estimate will always be normalized to sum to one
#' over the estimate grid.
#' @param x an sf or sfc object. The KDE grid will be estimated from the bounding box.
#' @param barrier An sf polygon object that defines an area where use is excluded.
#' @param norm Logical. Should each individual kernel be normalized to account for barriers
#' before addition to the total KDE. Defaults to \code{TRUE} if a barrier is specified.
#' @param ... Any arguments passed to \code{sf::st_make_grid} to create the KDE prediction grid. See \link[sf:st_make_grid]{sf::st_make_grid} for
#' description of additional arguments to make the KDE grid.
#' @import sf
#' @useDynLib crawlUtils, .registration = TRUE
#' @author Devin S. Johnson
#' @export
#'
cu_kde_ud <- function(x, barrier=NULL, norm, ...){
  grid <- sf::st_make_grid(x,...)
}

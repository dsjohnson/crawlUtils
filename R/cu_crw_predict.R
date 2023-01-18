#' @title Batch CRW Prediction for Multiple Animals
#' @description Uses a list of CRW fitted models and desired prediction times
#' to make location (and velocity) predictions for telemetered animals.
#' @param fit A CRW fit object
#' @param predTime A character string describing the desired frequency of prediction,
#' e.g., \code{predTime="1 hour"} or \code{predTime="15 min"}.
#' @param barrier An \code{sf} polygon object representing areas where the animal cannot access.
#' @param vis_graph A visibility graph constructed with the R package \code{pathroutr}, which is used
#' to reroute paths around barriers.
#' @param as_sf Logical. Return an \code{sf} points data frame (\code{TRUE}) or standard \code{crawl} prediction.
#' @param ... Additional arguments passed to the \code{\link[foreach]{foreach}} function, e.g.,
#' for error handling in the loop.
#' @details The R package \code{pathroutr} is necessary for use of the \code{barrier} rerouting.
#' it can be installed with the command
#' \code{install.packages('pathroutr', repos='https://jmlondon.r-universe.dev')}.
#' See 'https://github.com/jmlondon/pathroutr' for a description of use and constructing the
#' viability \code{vis_graph}.
#' @author Devin S. Johnson
#' @export
#' @import sf dplyr crawl foreach
#'
cu_crw_predict <- function(fit, predTime=NULL, barrier=NULL, vis_graph=NULL, as_sf=TRUE,...){
  locType <- NULL #handle 'no visible binding...'
  # progressr::handlers(global = TRUE)
  route <- !is.null(barrier) & !is.null(vis_graph)
  # p <- progressr::progressor(length(fit_list))
  pred <- crawl::crwPredict(fit, predTime=predTime, return.type="flat")
  if(route){
    if (!requireNamespace("pathroutr", quietly = TRUE)) stop("Please install {pathroutr} package: install.packages('pathroutr',repos='https://jmlondon.r-universe.dev')")
    pred <- pred %>% crawl::crw_as_sf(ftype="POINT", locType="p") %>% filter(locType=="p")
    pred <- pred %>% pathroutr::prt_trim(barrier)
    fix <- pathroutr::prt_reroute(pred, barrier, vis_graph, blend=FALSE)
    pred <- pathroutr::prt_update_points(fix, pred)
  }
  if(as_sf & !route) pred <- crw_as_sf(pred, ftype="POINT")
  # p()
  return(pred)
}

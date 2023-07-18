#' @title Batch CRW Posterior Path Simulation For Multiple Animals
#' @description Uses a list of CRW fitted models and desired simulation times
#' to make draws from the location (and velocity) posterior distribution for telemetered animals.
#' @param fit A CRW fit object
#' @param size The number of posterior draws. Defaults to 8 (See Details).
#' @param predTime A character string describing the desired frequency of prediction,
#' e.g., \code{predTime="1 hour"} or \code{predTime="15 min"}.
#' @param barrier An \code{sf} polygon object representing areas where the animal cannot access.
#' @param vis_graph A visibility graph constructed with the R package \code{pathroutr}, which is used
#' to reroute paths around barriers.
#' @param as_sf Logical. Return an \code{sf} points data frame list (\code{TRUE}) or standard \code{crawl} prediction list
#' @param ... Additional arguments passed to the \code{\link[foreach]{foreach}} function, e.g.,
#' for error handling in the loop.
#' @details The R package \code{pathroutr} is necessary for use of the \code{barrier} rerouting.
#' it can be installed with the command
#' \code{install.packages('pathroutr', repos='https://jmlondon.r-universe.dev')}.
#' See 'https://github.com/jmlondon/pathroutr' for a description of use and constructing the
#' viability \code{vis_graph}.
#' @author Devin S. Johnson
#' @export
#' @import sf dplyr foreach crawl
#'
cu_crw_sample <- function(fit, size=8, predTime=NULL, barrier=NULL, vis_graph=NULL, as_sf=TRUE,...){
  j <- fid <- NULL #handle 'no visible binding...'
  # progressr::handlers(global = TRUE)
  route <- !is.null(barrier) & !is.null(vis_graph)
  # p <- progressr::progressor(length(fit_list))
  simObj <- crawl::crwSimulator(fit, parIS = 0, predTime=predTime)
  out <- foreach(j=1:size)%do%{
    samp <- crawl::crwPostIS(simObj, fullPost = FALSE)
    attr(samp, "crw_type") <- "crwIS"
    if(route){
      if (! requireNamespace("pathroutr", quietly = TRUE)) stop("Please install {pathroutr}: install.packages('pathroutr',repos='https://jmlondon.r-universe.dev')")
      samp <- samp %>% crawl::crw_as_sf(ftype="POINT", locType="p")
      samp <- samp %>% pathroutr::prt_trim(barrier)
      fix <- pathroutr::prt_reroute(samp, barrier, vis_graph, blend=FALSE)
      samp$geometry[fix$fid] <- fix$geometry
      samp <- samp %>% dplyr::mutate(rep=j)
      attr(samp, "crw_type") <- "crwIS_sf"
    }
    if(as_sf & !route){
      samp <- crw_as_sf(samp,ftype="POINT") %>% dplyr::mutate(rep=j)
      attr(samp, "crw_type") <- "crwIS_sf"
    }
    if("fid" %in% colnames(samp)) samp <- dplyr::select(samp, -fid)
    samp
  }
  if(as_sf){
    attr(out, "crw_type") <- "crwIS_sf_list"
  } else{
    attr(out, "crw_type") <- "crwIS_list"
  }
  # p()
  return(out)
}




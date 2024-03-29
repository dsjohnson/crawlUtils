#' @title Calculate Effective Sample Size for a Set of CRW locations
#' @description Estimates the number of independent locations in a CRW data set
#' using the mutual information method of Bartoszek (2016).
#' @param fit A \code{crwFit} object (See \code{\link[crawl]{crwMLE}}).
#' @param aug Either a \code{\link[crawl]{crwPredict}} or \code{\link[crawl]{crwPostIS}} objects
#' from which the extra location times will be used in the calculation.
#' The \code{\link[crawl]{crw_as_sf}} transformed versions of these objects will also work.
#' @details This function uses the "mutual information" effective sample size of Bartoszek (2016) to calculate the
#' equivalent number of independent animal locations. It also calculates individual contributions of each location using the
#' regression effective sample size in Bartoszek (2016). The output is a named list with `Ne` equal to the overall sample
#' and `w` is a vector of weights that sum to 1 overall. If you want the ESS value of each observation `Ne * w` will provide it.
#' @references Bartoszek, K. (2016). Phylogenetic effective sample size.
#' Journal of Theoretical Biology. 407:371-386. (See https://arxiv.org/pdf/1507.07113.pdf).
#' @author Devin S. Johnson
#' @export
#'
cu_crw_ess <- function(fit, aug=NULL){
  df <- fit$data
  # v <- diag(var(df[,fit$coord]))
  # df <- fit$data[!duplicated(fit$data$TimeNum),]
  cf <- cu_crw_covfun(fit)
  if(!is.null(aug)){
    if(is.null(aug$TimeNum)) stop("The 'aug' argument does not appear to be a 'crwIS', 'crwPredict' or sf POINT geometry object.")
    times <- aug$TimeNum
    if(any(duplicated(times))) stop("There are duplicated times in the 'aug' object. Please remove these first.")
  } else{
    times <- df$TimeNum
    if(any(duplicated(times))) stop("There are duplicated times in the original data.")
  }
  n <- length(times)
  R <- cu_crw_covmat(times, corr=FALSE, cf=cf, E=times[1]-(times[n]-times[1]))

  Reig <- eigen(R, symmetric = TRUE)
  ln_det_V <- sum(log(Reig$values))
  sum_ln_det_Vj <- sum(log(diag(R)))
  ne <-  1 + (n-1)/log(exp(1) + (sum_ln_det_Vj - ln_det_V))
  w <- rep(NA, n)
  condR <- 1/diag(Reig$vectors %*% diag(1/Reig$values) %*% t(Reig$vectors))
  w <- (1/n) + ((n-1)/n)*(condR/diag(R))
  w <- w/sum(w)
  out <- list(Ne = ne, w=w)
  class(out) <- c("list","crwESS")
  return(out)
}

#' @title Calculate Effective Sample Size for a Set of CRW locations
#' @description Estimates the number of independent locations in a CRW data set
#' using the mutual information method of Bartoszek (2016).
#' @param fit A \code{crwFit} object (See \code{\link[crawl]{crwMLE}}).
#' @param aug Either a \code{\link[crawl]{crwPredict}} or \code{\link[crawl]{crwPostIS}} objects
#' from which the extra \code{predTime} location times will be used in the calculation.
#' The \code{\link[crawl]{crw_as_sf}} transformed versions of these objects will also work.
#' @details Uses the "mutual information" formulation of Bartoszek (2016) to calculate the
#' equivalent number of independent animal locations.
#' @references Bartoszek, K. (2016). Phylogenetic effective sample size.
#' Journal of Theoretical Biology. 407:371-386. (See https://arxiv.org/pdf/1507.07113.pdf).
#' @author Devin S. Johnson
#' @export
#'
cu_crw_ess <- function(fit, aug=NULL){
  df <- fit$data
  v <- diag(var(df[,fit$coord]))
  df <- fit$data[!duplicated(fit$data$TimeNum),]
  if(is.null(aug)){
    cf <- cu_crw_covfun(fit)
    times <- df$TimeNum
    n <- length(times)
    R <- cu_crw_covmat(times[-1], corr=FALSE, cf=cf, E=times[1])
    Sx <- matrix(v[1],n,n)
    Sx[2:n,2:n] <- Sx[2:n,2:n] + R
    Sy <- matrix(v[2],n,n)
    Sy[2:n,2:n] <- Sy[2:n,2:n] + R
    # Sx <- cov2cor(Sx)
    # Sy <- cov2cor(Sy)
  } else{
    if(!"TimeNum"%in%colnames(aug)) stop("The 'aug' argument is not the correct class. See ?cu_crw_ess.")
    cf <- cu_crw_covfun(fit)
    first_obs <- df$TimeNum[1]==aug$TimeNum[1]
    if(!first_obs){
      times <- c(df$TimeNum[1],aug$TimeNum)
    } else{
      times <- aug$TimeNum
    }
    times <- times[!duplicated(times)]
    n <- length(times)
    R <- cu_crw_covmat(times[-1], corr=FALSE, cf=cf, E=times[1])
    Sx <- matrix(v[1],n,n)
    Sx[2:n,2:n] <- Sx[2:n,2:n] + R
    Sy <- matrix(v[2],n,n)
    Sy[2:n,2:n] <- Sy[2:n,2:n] + R
    if(!first_obs){
      Sx <- Sx[-1,-1]
      Sy <- Sy[-1,-1]
      n <- n-1
    }
    # Sx <- cov2cor(Sx)
    # Sy <- cov2cor(Sy)
  }
  # ess <- (sum(solve(Sx,rep(1,nrow(Sx)))) + sum(solve(Sy,rep(1,nrow(Sy)))))/2
  ln_det_V <- sum(log(eigen(Sx)$values)) + sum(log(eigen(Sy)$values))
  sum_ln_det_Vj <- 0
  for(j in 1:n){
    sum_ln_det_Vj <- sum_ln_det_Vj + (log(Sx[j,j]) + log(Sy[j,j]))
  }
  ess <-  1 + (n-1)/log(exp(1) + 0.5*(sum_ln_det_Vj - ln_det_V))
  return(ess)
}

#' Estimate Confidence Interval Using Likelihood Ratio
#'
#' Estimate the \code{(1-a)}-level confidence interval of each parameter's estimate
#' from class \code{siestim}.
#'
#' @details
#' For an estimate \eqn{\theta}, the interval is estimated such that \cr
#' \deqn{}
#' where \eqn{z_a} is \code{(1-a)}-th quantile of the standard Normal distribution,
#' and \eqn{se(\theta)} is the estimated standard error of \eqn{\theta}.
#'
#' @param x object of class \code{siestim}.
#' @param level confidence interval level; default is 0.95.
#'
#' @return confidence interval matrix.

#'
#' @examples
getci2.siestim <- function(x, level = .95){
  #confidence interval using likelihood ratio
  chi2 <- qchisq(level, 4)
  intv <- x$logll - chi2

  fn <- function(par){
    # put constraints (important when doing optimization and computing hessian)
    if(par[1]<=0) par[1] <- 1e-5 # mean> 0
    if(par[2]<=0) par[2] <- 1e-5 # sd>0
    if(par[3]<=0) par[3] <- 1e-5 # pi>0; dgeom is NaN when pi<=0
    if(par[4]<0) par[4] <- 0 # w>=0
    if(par[3]>1) par[3] <- 1 #pi<=1; dgeom is NaN when pi>1
    if(par[4]>1) par[4] <- 1 # 0<=w<=1

    (-logll(par, unlist(x$data), x$control) - intv)^2
  }

  # find lower & upper bound
  low <- dfoptim::nmkb(x$par/2, fn, rep(0,4), x$par)
  upp <- dfoptim::nmkb(c(x$par[1:2]*3/2, (x$par[3:4]+1)/2), fn, x$par, c(Inf, Inf, 1, 1))

  ci <- cbind(low$par, upp$par)
  ci[,1] <- ifelse(ci[,1]<0, 0, ci[,1]) #if the lower bound is negative, make it zero
  ci[3:4,2] <- ifelse(ci[3:4,2]>1, 1, ci[3:4,2]) #if the upper bound of pi/w exceed 1, make it 1

  colnames(ci) <- c("lower", "upper")
  row.names(ci) <- names(x$par)

  return(ci)
}

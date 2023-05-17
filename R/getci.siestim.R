#' Estimate Confidence Interval Using Wald Interval
#'
#' Estimate the \code{(1-a)}-level confidence interval of each parameter's estimate
#' from class \code{siestim}.
#'
#' @details
#' For an estimate \eqn{\theta}, the interval is estimated as follows \cr
#' \deqn{\theta - z_{a/2}*se(\theta), \theta - z_{a/2}*se(\theta)}
#' where \eqn{z_a} is \code{(1-a)}-th quantile of the standard Normal distribution,
#' and \eqn{se(\theta)} is the estimated standard error of \eqn{\theta}.
#'
#' @param x object of class \code{siestim}.
#' @param level confidence interval level; default is 0.95.
#'
#' @return confidence interval matrix.
#' @export
#'
#' @examples
getci2.siestim <- function(x, level = .95){
  #confidence interval using std normal dist.
  z <- qnorm(level + (1 - level)/2)
  se <- x$se
  ci <- cbind(x$par-z*se, x$par+z*se)
  ci[,1] <- ifelse(ci[,1]<0, 0, ci[,1]) #if the lower bound is negative, make it zero
  ci[3:4,2] <- ifelse(ci[3:4,2]>1, 1, ci[3:4,2]) #if the upper bound of pi/w exceed 1, make it 1

  colnames(ci) <- c("lower", "upper")
  row.names(ci) <- names(x$par)

  return(ci)
}

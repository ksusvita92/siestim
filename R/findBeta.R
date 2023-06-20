#' Find Beta-distributed Parameters
#'
#' @description Search for Beta-distributed parameters given quantile and probability.
#'
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param log.p logical; if TRUE, \code{p} is returned as log(p).
#'
#' @return Return the Beta-distributed parameters, shape1 and shape2.
#' @export
#'
#' @examples
#' par <- c(10, 46)
#' p <- c(.3, .5, .9)
#' q <- qbeta(p, par[1], par[2])
#' findBeta(q, p)
#'
findBeta <- function(q, p, log.p = FALSE){
  # check
  if(length(q) != length(p)) stop("q and p have different length.")
  if(length(q) < 2) stop("Provide at least 2 quantiles.")

  # Define the objective function
  objFn <- function(par, q, p, log.p){
    phat <- pbeta(q, par[1], par[2], log.p = log.p)
    sum((phat - p)^2)
  }


  # find solution
  nlminb(rep(1,2), objFn, q = q, p = p, log.p = log.p, lower = 1e-300, upper = Inf)$par
}

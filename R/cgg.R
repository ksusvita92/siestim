#' @rdname cgg
#' @aliases rcgg, dcgg, pcgg
#' @title Compound Geometric Gamma
#'
#' @description Random generator, density, and cumulative density of a Compound Geometric
#' Gamma (CGG) distribution.
#'
#'
#' @param n number of observations.
#' @param x vector of quantiles representing the period of symptom onset times between
#' a pair of infector-infectee.
#' @param mu mean parameter of the serial interval distribution.
#' @param sigma standard deviation parameter of the serial interval distribution.
#' @param pi success probability parameter to sample the secondary case. \eqn(0 < pi <=1).
#'
#'
#' @details The distribution is to model the observed serial interval distribution coming from
#' primary-secondary transmission or indirect primary-secondary transmission
#' with some unknown intermediate unsampled individuals in between.
#'
#'
#' @return
#' @export
#'
#' @examples
#'
rcgg <- function(n,mu,sigma,pi){
  # use shape and rate
  a <- (mu/sigma)^2; b <- mu/sigma^2
  m <- rgeom(n,pi)
  x <- rgamma(n, (m+1)*a, b)
  return(x)
}
dcgg <- function(x, mu, sigma, pi){
  gt <- function(xi,mu,sigma,pi){
    # use shape and rate
    a <- (mu/sigma)^2; b <- mu/sigma^2

    #Find maximum m within tol=1e-10
    k <- 0; tol <- 1e-10
    y <- 1/gamma(a)
    v <- b^a * xi^a * (1-pi)
    while(y>tol && y<Inf){
      k <- k+1
      y <- tryCatch(v^k/gamma((k+1)*a), error=function(e) return(0), warning=function(w) return(0))
      if(is.nan(y)) y <- 0 #break the loop
      #NaN happens when gamma((m+1)*a)=Inf or when it's not defined
    }
    if(k>1) mmax <- k-1
    else mmax <- 0
    f <- sapply(0:mmax, function(i) dgamma(xi, (i+1)*a, b) * dgeom(i, pi))
    return(sum(f))
  }
  ft <- sapply(1:length(x), function(i) gt(xi=x[i], mu, sigma, pi))
  ft[is.infinite(ft)] <- 0 # it evaluate x=0 as Inf.
  return(ft)
}
pcgg <- function(x, mu, sigma, pi){
  pt <- function(xi,mu,sigma,pi){
    int <- integrate(f=dcgg, lower=0, upper=xi, mu=mu, sigma=sigma, pi=pi)
    val <- int$value
    return(val)
  }
  Ft <- sapply(1:length(x), function(i) pt(x[i], mu, sigma, pi))
  Ft[which(Ft>1)] <- 1
  return(Ft)
}

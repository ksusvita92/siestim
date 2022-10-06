#' Compound Geometric Gamma
#'
#' @description Random generator, density, cumulative density, and quantile functions
#' of a Compound Geometric
#' Gamma (CGG) distribution.
#'
#' @details The distribution is to model the observed serial interval distribution coming from
#' non-coprimary transmission (a primary-secondary transmission or indirect primary-secondary
#' transmission
#' with some unknown intermediate unsampled individuals in between). \cr\cr
#' The density function is given as follows \cr
#' \deqn{f(t) = \sum_m g(t|m) * Pr(M=m), t > 0, m = 0,1,...}
#' \eqn{g(t|m)} is the Gamma distributed density
#' function with shape \eqn{(m+1)\alpha} and rate \eqn{\beta}, and \eqn{Pr(M = m)} is the
#' probability mass function of Geometric distribution. Here, \eqn{\alpha, \beta} are shape
#' and rate parameters of the true serial interval distribution. \cr\cr
#' The quantile function is determined by solving the optimization function as follows
#' \deqn{arg min_t (F(t)-p)^2,}
#' where \eqn{F(t)} is the cumulative density function of CGG distribution and \eqn{p}
#' is a given probability.
#'
#' @param n number of observations.
#' @param x vector of quantiles representing the period of symptom onset times between
#' a pair of infector-infectee.
#' @param p vector of probabilities.
#' @param mu mean of the serial interval distribution.
#' @param sigma standard deviation of the serial interval distribution.
#' @param pi success probability to sample the secondary cases.
#' The value must be in (0,1].
#' @param log,log.p logical; if \code{TRUE}, the probabilities/densities \eqn{p} are returned
#' as \eqn{log(p)}.
#' @param lower.tail logical; if \code{TRUE}, the probabilities are returned as \eqn{Pr(X \le x)},
#' otherwise \eqn{Pr(X > x)}.
#'
#' @return \code{rcgg} generates random variables, \code{dcgg} returns the density,
#' \code{pcgg} returns the cumulative density, and \code{qcgg} returns
#' the quantile given probability \code{p} of CGG distribution. \cr\cr
#' The length of the result is determined by \code{n} for \code{rcgg}, the length of
#' \code{x} for \code{dcgg} and \code{pcgg}, and the length of \code{p} for \code{qcgg}.
#'
#' @examples x <- rcgg(10, 4.5, 2, .8)
#' fx <- dcgg(x, 4.5, 2, .7, T)
#' Fx <- pcgg(x, 4.5, 2, .7)
#'
#' @rdname cgg
#' @aliases rcgg
#' @export
rcgg <- function(n,mu,sigma,pi){
  # use shape and rate
  a <- (mu/sigma)^2; b <- mu/sigma^2
  m <- rgeom(n,pi)
  x <- rgamma(n, (m+1)*a, b)
  return(x)
}

#' @rdname cgg
#' @aliases dcgg
#' @export
dcgg <- function(x, mu, sigma, pi, log = FALSE){
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

  # compute the density
  ft <- sapply(1:length(x), function(i) gt(xi=x[i], mu, sigma, pi))
  ft[is.infinite(ft)] <- 0 # it evaluate x=0 as Inf.

  # is it log?
  if(log) ft <- log(ft)

  return(ft)
}

#' @rdname cgg
#' @aliases pcgg
#' @export
pcgg <- function(x, mu, sigma, pi, lower.tail = TRUE, log.p = FALSE){
  pt <- function(xi,mu,sigma,pi){
    int <- integrate(f=dcgg, lower=0, upper=xi, mu=mu, sigma=sigma, pi=pi)
    val <- int$value
    return(val)
  }

  # the cdf
  Ft <- sapply(1:length(x), function(i) pt(x[i], mu, sigma, pi))
  Ft[which(Ft>1)] <- 1

  # not lower tail?
  if(!lower.tail) Ft <- 1-Ft

  # log?
  if(log.p) Ft <- log(Ft)

  return(Ft)
}



#' @rdname cgg
#' @aliases qcgg
#' @export
qcgg <- function(p, mu, sigma, pi, lower.tail = TRUE, log.p = FALSE){
  obj <- function(x, quantile, mu, sigma, pi, lower.tail, log.p){
    (pcgg(x, mu, sigma, pi, lower.tail, log.p) - quantile)^2
  }

  res <- sapply(p, function(i) nlminb(mu, obj, quantile = i, mu = mu, sigma = sigma, pi = pi, lower.tail = lower.tail, log.p = log.p, lower = 1e-300, upper = Inf)$par)

  return(res)
}

#' Folded Gamma Difference
#'
#' @description Random generator, density, cumulative density, and quantile functions of
#' a Folded Gamma Difference (FGD) distribution.
#'
#' @details The distribution is to model the observed serial interval distribution coming from
#' primary-coprimary transmission, i.e. a type of transmission where two cases are linked,
#' but both were infected by the unseen third. \cr\cr
#' The density comes from Folded Gamma Difference (FGD),
#' which is derived from Gamma Difference Distribution [1]. The function has non-negative
#' real number support and is expressed as follows \cr
#' \deqn{f(t) = 2 \int g(s)*g(s-t) ds,} \cr
#' for \eqn{t \ge 0} and \eqn{t \le s \le \infty}. Here, \eqn{g(t), t > 0} is the density of
#' the true serial interval following Gamma distribution having shape \eqn{\alpha} and rate
#' \eqn{\beta}. \cr\cr
#' The quantile function is determined by solving the optimization function as follows
#' \deqn{arg min_t (F(t)-p)^2,}
#' where \eqn{F(t)} is the cumulative density function of FGD distribution and \eqn{p}
#' is a given probability.
#'
#' @param n number of observations.
#' @param x vector of quantiles representing the period of symptom onset times between
#' a pair of infector-infectee.
#' @param p vector of probabilities.
#' @param mu mean parameter of the serial interval distribution.
#' @param sigma standard deviation parameter of the serial interval distribution.
#' @param log,log.p loggical; if \code{TRUE}, the probabilities/densities \eqn{p} are returned
#' as \eqn{log(p)}.
#' @param lower.tail logical; if \code{TRUE}, the probabilities are returned as \eqn{Pr(X \le x)},
#' otherwise \eqn{Pr(X > x)}.
#'
#' @return \code{rfgd} generates random variables, \code{dfgd} returns the density,
#' \code{pfgd} returns the cumulative density,  and \code{qfgd} returns the quantile
#' given probability \code{p} of FGD distribution. \cr\cr
#' The length of the result is determined by \code{n} for \code{rfgd}, the length of
#' \code{x} for \code{dfgd} and \code{pfgd}, and the length of \code{p} for \code{qfgd}.
#'
#' @examples x <- rfgd(10, 4.5, 2)
#' fx <- dfgd(x, 4.5, 2)
#' Fx <- pfgd(x, 4.5, 2)
#' q <- qfgd(Fx, 4.5, 2)
#'
#' @references
#' \enumerate{
#'    \item Bernhard Klar (2015) A note on gamma difference distributions,
#' Journal of Statistical Computation and Simulation, 85:18, 3708-3715,
#' DOI: 10.1080/00949655.2014.996566
#' }
#'
#' @rdname fgd
#' @aliases rfgd
#' @export
rfgd <- function(n, mu, sigma){
  # use shape and rate
  a <- (mu/sigma)^2; b <- mu/sigma^2

  x1 <- rgamma(n, a, b)
  x2 <- rgamma(n, a, b)
  return(abs(x1-x2))
}




#' @rdname fgd
#' @aliases dfgd
#' @export
dfgd <- function(x, mu, sigma, log = FALSE){
  gt <- function(ti, mu, sigma){
    a <- (mu/sigma)^2; b <- mu/sigma^2
    if(ti>=0){
      integrand <- function(u, t) dgamma(u,a,b)*dgamma(u-t,a,b)
      ht <- tryCatch(integrate(integrand, ti, Inf, t = ti)$value,
                     error = function(e) return(0))
    } else ht <- 0

    if(is.nan(ht)) ht <- 0 #happens because we may multiply 0 and Inf
    return(ht)
  }

  # the pdf
  ft <- sapply(1:length(x), function(i) gt(ti=x[i], mu, sigma))
  ft <- 2*ft

  # log?
  if(log) ft <- log(ft)

  return(ft)
}



#' @rdname fgd
#' @aliases pfgd
#' @export
pfgd <- function(x, mu, sigma, lower.tail = TRUE, log.p = FALSE){
  pt <- function(xi,mu,sigma){
    int <- integrate(f=dfgd, lower=0, upper=xi, mu=mu, sigma=sigma)
    val <- int$value
    return(val)
  }

  # the cdf
  Ft <- sapply(1:length(x), function(i) pt(x[i], mu, sigma))
  Ft[which(Ft>1)] <- 1

  # not lower tail?
  if(!lower.tail) Ft <- 1-Ft

  # log?
  if(log.p) Ft <- log(Ft)

  return(Ft)
}



#' @rdname fgd
#' @aliases qfgd
#' @export
qfgd <- function(p, mu, sigma, lower.tail = TRUE, log.p = FALSE){
  obj <- function(x, quantile, mu, sigma, lower.tail, log.p){
    (pfgd(x, mu, sigma, lower.tail, log.p) - quantile)^2
  }

  res <- sapply(p, function(i) nlminb(0, obj, quantile = i, mu = mu, sigma = sigma, lower.tail = lower.tail, log.p = log.p, lower = 0, upper = Inf)$par)

  return(res)
}

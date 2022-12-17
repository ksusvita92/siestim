#' The negative log-likelihood functions
#'
#' @param params
#' @param dt
#' @param ctr
#'
#' @return
logll <- function(params, dt, ctr = list()){
  # see if mle or map is used
  if(is.null(ctr$prior.pi)) fp <- 1 else fp <- ctr$prior.pi
  if(is.null(ctr$prior.w)) fw <- 1 else fw <- ctr$prior.w
  if(is.null(ctr$prior.mu)) fm <- 1 else fm <- ctr$prior.mu
  if(is.null(ctr$prior.sigma)) fs <- 1 else fs <- ctr$prior.sigma


  f1 <- dcgg(dt, params[1], params[2], params[3])
  f2 <- dfgd(dt, params[1], params[2])
  #avoid log(0), log(Inf); give some offsets
  f1[f1==0 | !is.finite(f1)] <- f2[f2==0 | !is.finite(f1)] <- 1e-300 #non-finite occurs when e.g dgamma(0,.1,2)
  ft <- params[4]*f1 + (1-params[4])*f2
  if(fm == 0 || !is.finite(fm)) fm <- 1e-300
  if(fs == 0 || !is.finite(fs)) fs <- 1e-300
  if(fp == 0 || !is.finite(fp)) fp <- 1e-300
  if(fw == 0 || !is.finite(fw)) fw <- 1e-300


  # result
  nll <- -sum(log(ft))-log(fp)-log(fw)-log(fm)-log(fs)

  return(nll)
}


#' The optimization function
#'
#' @param params0
#' @param dt
#' @param lower
#' @param upper
#' @param ctr
#'
#' @importFrom dfoptim nmkb
#' @importFrom numDeriv hessian
#'
#' @return
#' @export
opt <- function(params0, dt, lower, upper, ctr = list()){
  # omit NA
  dt <- dt[!is.na(dt)]


  # estimate parameters
  myestim <- dfoptim::nmkb(params0, logll, lower, upper, dt=dt)
  params <- myestim$par
  mylogll <- -myestim$value
  msg <- myestim$message


  # estimate var-cov matrix using hessian
  hessian <- numDeriv::hessian(logll, params, dt = dt, ctr = ctr)
  varcov <- solve(hessian) # the likelihood is already negative
  se <- sqrt(abs(diag(varcov)))

  names(params) <- names(se) <- c("mu", "sigma", "pi", "w")
  out <- list(par=params, se=se, logll=mylogll, msg=msg)
  return(out)
}


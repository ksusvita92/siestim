#' The negative log-likelihood functions
#'
#' @param params
#' @param dt
#' @param ctr
#' @param prior.pi
#' @param prior.w
#'
#' @return
logll <- function(params, dt, prior.pi = NULL, prior.w = NULL, ctr = list()){
  #setup
  #get error when computing gamma(0)
  if(params[1] == 0) params[1] <- 1e-300
  if(params[2] == 0) params[2] <- 1e-300
  if(params[3] == 0) params[3] <- 1e-300
  if(params[3] > 1) params[3] <- 1
  if(params[4] < 0) params[4] <- 0
  if(params[4] > 1) params[4] <- 1

  # see if mle or map is used
  if(is.null(prior.pi)) fp <- 1
  else fp <- dbeta(params[3], prior.pi[1], prior.pi[2])

  if(is.null(prior.w)) fw <- 1
  else fw <- dbeta(params[4], prior.w[1], prior.w[2])

  f1 <- dcgg(dt, params[1], params[2], params[3])
  f2 <- dfgd(dt, params[1], params[2])
  f1[f1==0 | is.infinite(f1) | is.na(f1)] <- 1e-300 #avoid log(0)
  f2[f2==0 | is.infinite(f2) | is.na(f2)] <- 1e-300 #avoid log(0)
  if(fp == 0) fp <- 1e-300
  if(fw == 0) fw <- 1e-300
  ft <- params[4]*f1 + (1-params[4])*f2

  # result
  nll <- -sum(log(ft)) - log(fp) - log(fw)

  return(nll)
}



#' The optimization function
#'
#' @param params0
#' @param dt
#' @param lower
#' @param upper
#' @param ctr
#' @param prior.pi
#' @param prior.w
#'
#' @importFrom dfoptim nmkb
#' @importFrom numDeriv hessian
#'
#' @return
#' @export
opt <- function(params0, dt, lower, upper, prior.pi = NULL, prior.w = NULL, ctr = list()){
  # omit NA
  dt <- dt[!is.na(dt)]

  # estimate parameters
  myestim <- dfoptim::nmkb(params0, logll, dt=dt, prior.pi = NULL, prior.w = NULL, control = ctr, lower=lower, upper=upper)

  params0 <- myestim$par
  mylogll <- -myestim$value
  convergence <- myestim$convergence

  # show message.
  # if estimates exceed bound(s), return with code -1
  if(all(params0[1:2] > lower[1:2]) && all(params0[1:2] < upper[1:2])){
    msg <- myestim$message
  } else{
    msg <- "Estimates exceed bound(s)."
    convergence <- -1
  }

  # estimate var-cov matrix using hessian
  hessian <- numDeriv::hessian(logll, params0, dt = dt, prior.pi = NULL, prior.w = NULL, ctr = ctr)
  varcov <- solve(hessian) # the likelihood is already negative

  param <- myestim$par
  names(param) <- c("mu", "sigma", "pi", "w")
  se <- sqrt(abs(diag(varcov)))
  names(se) <- c("mu", "sigma", "pi", "w")

  out <- list(par=param, se=se, varcov=varcov, logll=mylogll, convergence = convergence, msg=msg)
  return(out)
}

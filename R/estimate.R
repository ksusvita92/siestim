#' Serial Interval Estimation
#'
#' An optimization method (MLE/MAP) to estimate serial interval distribution.
#'
#' @details If the transmission networks are non-unique, \code{x} must be
#' expressed as a list, in which each element represents the observed serial interval data
#' per one transmission network. \cr\cr
#' \code{init} is a numeric vector (cannot be \code{NA} or \code{NULL}) which
#' (by order) depicts the initial value of \eqn{\mu, \sigma, \pi, w}. \eqn{\mu, \sigma} are
#' mean and standard deviation of the Gamma distributed serial interval, \eqn{\pi} is the
#' success probability to sample the secondary cases in non-coprimary transmissions, see
#' \code{?cgg}, and \eqn{w} is the probability of non-coprimary transmissions in the data. \cr\cr
#' \code{prior.pi}, \code{prior.w} supply the parameters of Beta distribution. By default,
#' \code{prior.pi}, \code{prior.w} are non-informative. If not \code{NULL}, MAP (maximum a
#' posteriori) method is used to estimate the parameters, otherwise MLE (maximum likelihood
#' estimation) method is used. \cr\cr
#' \code{control} argument is a list that can supply any of the following components:
#' \itemize{
#'    \item \code{tol} convergence tolerance; iteration is terminated when the absolute
#'    difference in function value between successive iteration is below \code{tol}.
#'    Default is 1.e-06.
#'    \item \code{regsimp} a logical variable indicating whether the starting parameter
#'    configuration is a regular simplex. Default is \code{TRUE}.
#'    \item \code{maximize} a logical variable indicating whether the objective function
#'    should be maximized. Default is \code{FALSE}.
#' } \cr
#' The observed serial intervals are assumed coming from two transmission types: coprimary and
#' non-coprimary transmission. The density is a mixture of the two transmission distribution,
#' which is given as follows
#' \deqn{f(t) = w*g(t|\mu,\sigma,\pi) + (1-w)*h(t|\mu,\sigma),  t \ge 0,}
#' where \eqn{g(t), h(t)} are the densities of Compound Geometric Gamma (CGG) distribution
#' and Folded Gamma Difference (FGD) distribution, respectively. The density above shares the
#' same parameters with the true serial interval distribution, i.e. the mean \eqn{\mu} and the
#' standard deviation \eqn{\sigma}.
#'
#' @param x vector/list of vectors of the observed symptom onset interval between a pair
#' of infector-infectee. See "Details".
#' @param init initial values of \eqn{\mu, \sigma, \pi, w}, where \eqn{\pi, \sigma > 0},
#' \eqn{0 < \pi \le 1} and
#' \eqn{0 \le w \le 1}. See "Details".
#' @param ci confidence interval level for each estimate. By default, it is 0.95.
#' @param lower,upper bounds on the parameters, cannot be \code{NA} or \code{NULL}. By order,
#' it is \eqn{\mu, \sigma, \pi, w}.
#' @param ncore number of physical CPUs/cores to be used during the computation. If \code{ncore > 1},
#' the evaluation is done by parallel.
#' @param prior.pi,prior.w vector of Beta distributed parameters for \eqn{\pi, w}. See "Details".
#' @param control a list of control parameters. See "Details".
#'
#' @import foreach
#' @import parallel
#' @import tcltk
#' @import doSNOW
#'
#' @return The function returns the following values: \cr
#' \itemize{
#'   \item \code{par}: point estimation of each parameter.
#'   \item \code{se}: standard error of each estimate.
#'   \item \code{ci}: the \code{ci x 100%} confidence interval of each estimate.
#'   \item \code{run.time}: total execution time.
#'   \item \code{success.run}: proportion of convergence estimations.
#'   \item \code{record}: a data frame that records the estimation results.
#' }
#'
#' @export
#'
#' @examples
#' data <- simulateData(100, 4, 2, R0 = 1)
#' ncore <- parallel::detectCores()
#' est <- estimate(x = data$transnet$onset_diff, init = c(2,2,.5,.5), lower = c(10,10,1,1), ncore = ncore)
#' est
estimate <- function(x, init, lower, upper, ci = .95, ncore = 1, prior.pi = NULL, prior.w = NULL, control = list()){
  # start the clock
  start <- Sys.time()


  # check bounds
  lower <- ifelse(lower < 0, 0, lower)
  upper[3:4] <- ifelse(upper[3:4] > 1, 1, upper[3:4])


  # check if x is in list format
  if(class(x) != "list"){
    if(class(x) %in% c("numeric", "double", "integer")){
      x <- list(x)
      warning("x is considered as a single data.", call. = F)
    } else stop("x is not of class 'numeric', 'double', or 'integer'.")
  }


  # check if ncore > 1
  if(ncore == 1){
    warning("Parallel is not performed. See 'Details' if you want otherwise.", call. = F)
  }


  #estimate one data
  tmpfn <- function(params0, dt, low, upp, prior.pi = NULL, prior.w = NULL, ctr){
    result <- opt(params0, dt, low, upp, prior.pi = NULL, prior.w = NULL, ctr)

    # get all outputs
    par <- result$par
    se <- result$se
    output <- data.frame(mu = par[1], sigma = par[2], pi = par[3], w = par[4],
                         se.mu = se[1], se.sigma = se[2], se.pi = se[3], se.w = se[4],
                         ll = result$logll, convergence = result$convergence, msg = result$msg, stringsAsFactors = FALSE)
    return(output)
  }


  # show progress bar
  pb <- tkProgressBar("SI estimation", "Progress...", 0, length(x), 0)
  progress <- function(n){
    info <- sprintf("%1.0f%% done", n/length(x)*100)
    setTkProgressBar(pb, n, "SI estimation", info)
  }
  opts <- list(progress=progress)


  # do parallel?
  if(ncore == 1){
    res <- foreach(i = 1:length(x),
                   .combine = rbind,
                   .export = c("dcgg", "dfgd", "logll", "opt"),
                   .options.snow=opts, .errorhandling="pass") %do% {
                     tryCatch(tmpfn(init, x[[i]], lower, upper, prior.pi = NULL, prior.w = NULL, control),
                              error = function(e){
                                msg <- conditionMessage(e)
                                output <- data.frame(mu = NA, sigma = NA, pi = NA, w = NA,
                                                     se.mu = NA, se.sigma = NA, se.pi = NA, se.w = NA,
                                                     ll = NA, convergence = NA, msg = msg, stringsAsFactors = FALSE)
                                return(output)
                              })
                   }
  } else{
    if(ncore > parallel::detectCores()) ncore <- parallel::detectCores()
    # setup parallel
    cl <- makeCluster(ncore, type = "SOCK")
    registerDoSNOW(cl)

    res <- foreach(i = 1:length(x),
                   .combine = rbind,
                   .export = c("dcgg", "dfgd", "logll", "opt"),
                   .options.snow=opts, .errorhandling="pass") %dopar% {
                     .GlobalEnv$dcgg <- dcgg
                     .GlobalEnv$dfgd <- dfgd
                     .GlobalEnv$logll <- logll
                     .GlobalEnv$opt <- opt
                     tryCatch(tmpfn(init, x[[i]], lower, upper, prior.pi = NULL, prior.w = NULL, control),
                              error = function(e){
                                msg <- conditionMessage(e)
                                output <- data.frame(mu = NA, sigma = NA, pi = NA, w = NA,
                                                     se.mu = NA, se.sigma = NA, se.pi = NA, se.w = NA,
                                                     ll = NA, convergence = NA, msg = msg, stringsAsFactors = FALSE)
                                return(output)
                              })
                   }

    stopCluster(cl)
  }
  close(pb)


  # get output
  tmp <- res
  res <- res[res$convergence == 0,]
  sc <- nrow(res)/length(x) # this is the rate of successfull runs
  if(nrow(res) == 0){
    res <- res
    warning("The estimates are not convergent or some are out of bounds. Try to increase the interval.", call. = F)
  }


  # get all the estimates
  par <- c(mean(res$mu), mean(res$sigma), mean(res$pi), mean(res$w))

  if(length(x) == 1){
    se <- c(res$se.mu, res$se.sigma, res$se.pi, res$se.w)
  } else{
    var <- c(var(res$mu)+mean(res$se.mu)^2, var(res$sigma)+mean(res$se.sigma)^2, var(res$pi)+mean(res$se.pi)^2, var(res$w)+mean(res$se.w)^2)
    se <- sqrt(var)
  }


  #confidence interval
  z <- qnorm(ci + (1 - ci)/2)
  cim <- cbind(par-z*se, par+z*se)
  cim[,1] <- ifelse(cim[,1]<0, 0, cim[,1]) # if the lower bound is negative, make it zero
  cim[3:4,2] <- ifelse(cim[3:4, 2] > 1, 1, cim[3:4,2]) # if the upper bound of pi/w exceed 1, make it 1



  # set param names
  names(par) <- names(se) <- rownames(cim) <- c("mu", "sigma", "pi", "w")
  colnames(cim) <- c("lower", "upper")



  output <- list(par = par, se = se, ci = cim, run.time = Sys.time()-start, success.run = sc, record = tmp)
  new_siestim(output)
}




#' @export
print.siestim <- function(x){
  cat("\n ----- Serial Interval Estimation ----- \n")
  cat("============================================")
  cat("\n \n")

  cat("Parameter estimations: \n")
  print(round(x$par, 3))
  cat("\n")

  cat("Standard error: \n")
  print(round(x$se, 3))
  cat("\n")

  cat("Confidence intervals: \n")
  print(round(t(x$ci), 3))
  cat("\n")

  cat("Running time: ", x$run.time, " \n")
  cat("Successful run: ", x$success.run*100, "% \n")
}

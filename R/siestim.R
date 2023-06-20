#' Serial Interval Estimation
#'
#' An optimization method (MLE/MAP) to estimate serial interval distribution.
#'
#' @details
#' If the transmission trees are non-unique, \code{x} must be
#' expressed as a list, in which each element represents the observed serial interval data
#' per one transmission tree. \cr
#'
#' If \code{x} is a list of length greater than one, the estimation will be performed by parallel.
#' By default, only 2 cores will be used.
#'
#' \code{init} is a numeric vector (cannot be \code{NA} or \code{NULL}) which
#' (by order) depicts the initial value of \code{mu}, \code{sigma}, \code{pi}, \code{w}.
#' \itemize{
#'    \item \code{mu}, \code{sigma} are mean and standard deviation of the Gamma distributed serial interval,
#'    \item \code{pi} is the success probability to sample the secondary cases in non-coprimary transmissions (see \code{?cgg}),
#'    \item \code{w} is the probability of non-coprimary transmissions in the data.
#' }
#' \cr
#' \code{control} argument is a list that can supply any of the following components:
#' \itemize{
#'    \item \code{prior.mu}: a function if the prior distribution of \code{mu} is provided.
#'    \item \code{prior.sigma}: a function if the prior distribution of \code{sigma} is provided.
#'    \item \code{prior.pi}: a function if the prior distribution of \code{pi} is provided.
#'    \item \code{prior.w}: a function if the prior distribution of \code{w} is provided.
#'    \item \code{ncore}: number of physical CPUs/cores to be used during the computation.
#'    Only when \code{length(x)} > 1. By default \code{ncore} = 2.
#' } \cr
#' The observed serial intervals are assumed coming from two transmission types: coprimary and
#' non-coprimary transmissions. The density is a mixture of the two transmission distributions,
#' which is given as follows
#' \deqn{f(t) = w*g(t|mu,sigma,pi) + (1-w)*h(t|mu,sigma),  t \ge 0,}
#' where \eqn{g(t), h(t)} are the densities of Compound Geometric Gamma (CGG) distribution
#' and Folded Gamma Difference (FGD) distribution, respectively. The density above shares the
#' same parameters with the true serial interval distribution, i.e. the mean \code{mu} and the
#' standard deviation \code{sigma}. \cr\cr
#' If the transmission trees are not unique (i.e. \code{x} is a list of length > 1), the estimations
#' are computed by taking average over all estimates from the transmission trees. \cr
#'
#' @param x vector or a list of vectors of the observed symptom onset interval between a pair
#' of infector-infectee. See 'Details'.
#' @param init initial values of the estimates (\code{mu}, \code{sigma}, \code{pi}, \code{w}), where \code{mu}, \code{sigma} > 0,
#' 0 < \code{pi} <= 1 and
#' 0 <= \code{w} <= 1. See 'Details'.
#' @param lower,upper bounds on the estimates, cannot be \code{NA} or \code{NULL}. By order,
#' it is \code{mu}, \code{sigma}, \code{pi}, \code{w}.
#' @param control a list of control parameters. See 'Details'.
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
#'   \item \code{logll}: log-likelihood value evaluated at \code{par}.
#'   \item \code{data}: serial interval data.
#'   \item \code{msg}: convergence message of the optimization.
#'   \item \code{record}: data frame representation of the estimates and their standard errors.
#' }
#'
#' @export
#'
#' @examples
#' #simulate data
#' set.seed(1)
#' n <- 500
#' par <- c(5, 2, .4, .7)
#' m <- rbinom(1, n, par[4])
#' data <- lapply(1:5, function (i) c(rcgg(m,par[1],par[2],par[3]), rfgd(n-m,par[1],par[2])))
#'
#' #estimate the parameter
#' est <- siestim(data, c(7, 3, .5, .5), rep(0, 4), c(10, 5, 1, 1))
#' est
#' plot(est)
#'
siestim <- function(x, init, lower, upper, control = list(ncore = 2)){
  # check bounds
  lower <- ifelse(lower < 0, 0, lower)
  upper[3:4] <- ifelse(upper[3:4] > 1, 1, upper[3:4])


  if(is.numeric(x) || (is.list(x) && length(x)==1)){
    nx <- 1
    # run the optimization
    if(is.list(x)) x <- x[[1]]
    output <- tryCatch(opt(init, x, lower, upper, control),
                    error = function(e){
                      par <- se <- rep(NA,4)
                      names(par) <- names(se) <- c("mu", "sigma", "pi", "w")
                      return(list(par=par, se=se, logll = NA, msg=conditionMessage(e)))
                    })

    output$control <- control
    output$data <- x
    output$record <- cbind(as.data.frame(as.list(c(output$par, output$se))), logll = output$logll, msg = output$msg)
    names(output$record)[5:8] <- paste("se", names(output$se), sep = ".")
  } else if(is.list(x) && length(x)>1){
    nx <- length(x)
    #create text progress bar
    pb <- txtProgressBar(max = length(x), title = "siestim", label = "Progress...", style = 3)
    prog <- function(n) setTxtProgressBar(pb, n)

    #set up clusters
    cl <- makeCluster(control$ncore, type = "SOCK")
    registerDoSNOW(cl)

    res <- foreach(i=1:length(x), .combine = rbind,
                   .export = c("dcgg", "dfgd", "logll", "opt"),
                   .errorhandling="pass", .options.snow = list(progress=prog)) %dopar% {
                     .GlobalEnv$dcgg <- dcgg
                     .GlobalEnv$dfgd <- dfgd
                     .GlobalEnv$logll <- logll
                     .GlobalEnv$opt <- opt

                     myest <- tryCatch(opt(init, x[[i]], lower, upper, control),
                                       error = function(e){
                                         par <- se <- rep(NA,4)
                                         names(par) <- names(se) <- c("mu", "sigma", "pi", "w")
                                         return(list(par=par, se=se, logll=NA, msg=conditionMessage(e)))
                                       })
                     names(myest$se) <- paste("se", names(myest$se), sep = ".")
                     out <- as.data.frame(as.list(c(myest$par, myest$se)))
                     out$logll <- myest$logll
                     out$msg <- myest$msg
                     out
                   }

    stopCluster(cl)


    # get convergence estimates and their variance
    #tmp <- res[res$msg == "Successful convergence",]
    tmp <- res
    par <- sapply(tmp[,1:4], function(x) mean(x, na.rm = T))
    se <- sqrt(sapply(tmp[,5:8]^2, function(x) mean(x, na.rm = T)) + sapply(tmp[,1:4], function(x) var(x, na.rm = T)))
    names(se) <- names(par)
    loglike <- logll(par, unlist(x), control)/nx

    output <- list(par = par, se = se, logll = -loglike, control = control, msg = res$msg, record = res)
  }

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

  if(nrow(x$record)==1){
    cat("Convergence message: \n")
    print(x$msg)
    cat("\n")
  }

}

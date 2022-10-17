#' Simulate Transmission Data
#'
#' Simulate transmission data.
#' The serial interval between two connected cases are generated
#' from two distributions: Compound Geometric Gamma (CGG) and Folded Gamma Difference (FGD).
#'
#' @details The underlying transmissions in the data are classified as
#' non-coprimary and coprimary transmissions. The former refers to a direct/indirect
#' transmission between primary case to his/her sampled secondary case. The latter refers to
#' a transmission in which both primary and secondary cases were infected by an unseen
#' (unsampled) case. \code{w} represents the fraction of non-coprimary transmission
#' in the transmsission network. \cr\cr
#'
#' @param n number of data to be generated.
#' @param mu mean of the serial interval.
#' @param sigma standard deviation of the serial interval.
#' @param pi success probability to sample secondary cases. \code{0<pi<=1}.
#' @param w probability of non-coprimary transmissions. See "Details".
#' @param kappa expected number of contacts per case. If \code{pi > 0}, the expected number
#' of contacts is modelled by Poisson distribution with mean \code{kappa*pi}.
#'
#' @return A list of data frames which include the epidata and the transmission network.
#' The epidata includes the following information:
#' \itemize{
#'   \item \code{inf.ID,inf.source} case IDs and their sampled infectors. If \code{pi>0},
#'   the transmission may consist of few intermediate unsampled cases.
#'   \item \code{inf.times} time when a case developed symptoms.
#' } \cr
#' The transmission network includes the following information:
#' \itemize{
#'   \item \code{case_i,case_j} infector and infectee IDs.
#'   \item \code{ti,tj} symptom onset times of the infector and infectee.
#'   \item \code{onset_diff} observed serial interval, computed by \code{tj-ti}. Assumed to
#'   be non-negative.
#' }
#'
#' @export
#'
#' @examples
#' data <- simulateData(100, 4, 2, kappa = 1)
#' plot(data)
simulateData <- function(n, mu, sigma, pi = .5, w = .7, kappa = 3){
  # generate observed serial interval with par
  n1 <- rbinom(1, n, w)
  x1 <- rcgg(n1, mu, sigma, pi)
  x2 <- rfgd((n-n1), mu, sigma)
  si <- c(x1,x2)



  # generate contacts
  case_i <- case_j <- c()
  pools <- c()
  ti <- tj <- c()
  nI <- 0 # acrually at t=0, we have ID_0, this is for indexing only
  nSus <- n
  #



  # generate and sample descendant of ID_0
  a <- max(round(rpois(1, kappa) * pi), 1)
  case_i <- c(case_i, rep("ID_0", a))
  case_j <- c(case_j, paste("ID", (nI+1):(nI+a), sep = "_"))
  ti <- rep(0, a)
  tj <- si[1:a]
  pools <- case_j
  #



  # update dynamic
  nI <- nI + a
  nSus <- nSus - a
  #



  while(nSus > 0){
    if(length(pools) == 0){
      nI <- nI + 1
      inf <- paste("ID", nI, sep = "_") # import case
      a <- max(round(rpois(1, kappa) * pi), 1)
    } else{
      inf <- sample(pools, 1)
      a <- round(rpois(1, kappa) * pi)
    }

    # update case i and j
    pools <- pools[pools != inf]

    if(nSus < a) a <- nSus # to make sure we can only sample nSample cases

    #generate & sample decendent of inf
    if(a > 0){
      case_j <- c(case_j, paste("ID", (nI+1):(nI+a), sep = "_"))
      case_i <- c(case_i, rep(inf, a))

      if(inf %in% case_j) tmp <- tj[which(case_j == inf)]
      else tmp <- abs(max(tj)-abs(rnorm(1)))
      ti <- c(ti, rep(tmp, a))
      tj <- c(tj, tmp + si[which(case_i == inf)])

      pools <- c(pools, paste("ID", (nI+1):(nI+a), sep = "_"))
      nI <- nI + a
      nSus <- nSus - a
    }
  }
  tc <- data.frame(case_i = case_i, case_j = case_j, ti = ti, tj = tj, onset_diff = si)
  epi <- data.frame(inf.ID = case_j, inf.source = case_i, inf.times = tj)
  #



  # output
  out <- list(epidata = epi, transnet = tc)
  new_simdata(out)
}



#' @export
print.simdata <- function(x){
  cat("\n --- Showing the first six generated cases --- \n")
  cat("============================================== \n")
  print(x$epidata[1:6,])
  cat("\n")
}




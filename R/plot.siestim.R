#' Plot Serial Interval Estimation
#'
#' @description Draw plots of serial interval distribution, observed serial interval
#' distribution, and each parameter estimate.
#'
#' @return The function returns the following plots:
#' \itemize{
#'    \item \code{Serial Interval Distribution} depicts the distribution of the Gamma-distributed
#'    serial interval with mean \code{mu} and standard deviation \code{sigma}. If the transmission
#'    networks are non-unique, the distribution
#'    of each network will be drawn in gray color, and the final estimated distribution will be drawn in
#'    coral color.
#'    \item \code{Mixture Distribution} depicts the densities of the mixture model and
#'    its mixture components (coprimary and non-coprimary).
#'    \item \code{Distributions of Estimates} shows figure of the distribution of each
#'    estimate. If the transmission networks are non-unique, a histogram will be drawn.
#' }
#'
#' @param x object of class \code{siestim}
#'
#' @import ggplot2
#' @import patchwork
#'
#' @export
plot.siestim <- function(x){
  rec <- x$record
  rec <- na.omit(rec)
  par <- x$par


  # make serial interval plot
  alpha <- (rec$mu/rec$sigma)^2
  beta <- rec$mu/rec$sigma^2

  # draw empty plot
  plsi <- ggplot() + theme_light()

  # draw all dist in record
  for(i in 1:nrow(rec)){
    maxsi <- qgamma(.99, alpha[i], beta[i])
    minsi <- qgamma(.01, alpha[i], beta[i])

    maxsi <- maxsi+2
    minsi <- max(0, minsi-2)

    dt <- data.frame(x = seq(minsi, maxsi, length.out = 100))
    dt$y <- dgamma(dt$x, alpha[i], beta[i])
    plsi <- plsi + geom_line(aes(x = x, y = y, col = "est. per data"), dt)
  }

  # draw the estimate
  maxsi <- qgamma(.99, (par[1]/par[2])^2, par[1]/par[2]^2)
  minsi <- qgamma(.01, (par[1]/par[2])^2, par[1]/par[2]^2)

  maxsi <- maxsi+2
  minsi <- max(0, minsi-2)

  dt <- data.frame(x = seq(minsi, maxsi, length.out = 100))
  dt$y <- dgamma(dt$x, (par[1]/par[2])^2, par[1]/par[2]^2)
  plsi <- plsi + geom_line(aes(x = x, y = y, col = "final est."), dt) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 9)) +
    scale_color_manual(values = c("est. per data" = "gray", "final est." = "coral"), name = "") +
    labs(x = "", y = "density", title = "Gamma-distributed Serial Interval")
  #



  # make observed serial interval plot
  cggsi <- qcgg(.99, par[1], par[2], par[3])
  fgdsi <- qfgd(.99, par[1], par[2])
  cggsi2 <- qcgg(.01, par[1], par[2], par[3])
  fgdsi2 <- qfgd(.01, par[1], par[2])

  maxmix <- max(cggsi, fgdsi) + 2
  minmix <- max(0, min(cggsi2, fgdsi2)-2)

  dt <- data.frame(x = seq(minmix, maxmix, length.out = 100))
  dt$cgg <- dcgg(dt$x, par[1], par[2], par[3])
  dt$fgd <- dfgd(dt$x, par[1], par[2])

  plosi <- ggplot(dt) +
    geom_line(aes(x = x, y = par[4]*cgg, linetype = "non-cop. trans.")) +
    geom_line(aes(x = x, y = (1-par[4])*fgd, linetype = "cop. trans.")) +
    geom_line(aes(x = x, y = par[4]*cgg + (1-par[4])*fgd), col = "coral") +
    theme_light() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 10)) +
    scale_linetype_manual(values = c("non-cop. trans." = "dashed", "cop. trans." = "dotdash"), name = "") +
    labs(x = "", y = "", title = "Mixture Density")
  #



  # plot the estimate record
  #if(is.null(x$prior.pi)) p.pi <- rep(0, 100)
  #else p.pi <- dbeta(seq(0,1,length.out = 100), x$prior.pi[1], x$prior.pi[2])

  #if(is.null(x$prior.w)) p.w <- rep(0,100)
  #else p.w <- dbeta(seq(0,1,length.out = 100), x$prior.w[1], x$prior.w[2])

  plmu <- ggplot(rec) +
    geom_histogram(aes(x = mu, y = ..density..), col="black", fill="white", bins = 30) +
    geom_vline(aes(xintercept = par[1]), col = "coral", lty = "dashed") +
    theme_light() +
    labs(x = "days", y = expression(hat(mu)))
  plsg <- ggplot(rec) +
    geom_histogram(aes(x = sigma, y = ..density..), col="black", fill="white", bins = 30) +
    geom_vline(aes(xintercept = par[2]), col = "coral", lty = "dashed") +
    theme_light() +
    labs(x = "days", y = expression(hat(sigma)), title = "")
  plpi <- ggplot() +
    geom_histogram(aes(x = pi, y = ..density..), rec, col="black", fill="white", bins = 30) +
    geom_vline(aes(xintercept = par[3]), col = "coral", lty = "dashed") +
    #geom_line(aes(x = seq(0,1,length.out=100), y = p.pi)) +
    theme_light() +
    labs(x = "", y = expression(hat(pi)), title = "")
  plw <- ggplot() +
    geom_histogram(aes(x = w, y = ..density..), rec, col="black", fill="white", bins = 30) +
    geom_vline(aes(xintercept = par[4]), col = "coral", lty = "dashed") +
    #geom_line(aes(x = seq(0,1,length.out=100), y = p.w)) +
    theme_light() +
    labs(x = "", y = expression(hat(w)), title = "")


  # output
  out <- list((plsi+plosi), (plmu+plsg)/(plpi+plw))
  new_plotsi(out)
}

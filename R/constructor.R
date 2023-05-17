#' Define class "siestim"
new_siestim <- function(x) structure(x, class = "siestim")



#' Define class plotsi
new_plotsi <- function(x) structure(x, class = "plotsi")



#' Create generic function getci
#' @export
getci2 <- function(x, level = .95) UseMethod("getci2", x)
getci <- function(x, level = .95) UseMethod("getci", x)

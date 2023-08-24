#' Define class "siestim"
new_siestim <- function(x) structure(x, class = "siestim")



#' Define class plotsi
new_plotsi <- function(x) structure(x, class = "plotsi")



#' Define class outbreak
new_simOutbreak <- function(x) structure(x, class = "simOutbreak")


#' Create generic function getci
#' @export
getci <- function(x, level = .95) UseMethod("getci", x)
getci2 <- function(x, level = .95) UseMethod("getci2", x)


#' Create generic function createTC
#' @export
createTC <- function(x, downsampling_prob, cutoff, dna_model = "N", seed = 1) UseMethod("createTC", x)

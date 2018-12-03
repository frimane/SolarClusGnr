
#' Extraterrestrial horizontal solar irradiance generation
#'
#' Generic function to generate the extraterrestrial horizontal solar irradiance data.
#' It can genarate at any time resolution.
#'
#' @usage  extEHIcalc(phi = 0, lg = 0, tStep = 60)
#'
#' @param phi numreic(1) between -90 and 90 represents the latitude of the needed location.
#' @param lg numreic(1) between -180 and 180 represents the longitude of the needed location.
#' @param tStep numreic(1) represents the needed time step in seconds. Example to generate 1-min
#' data, tStep = 60.
#'
#' @return a matrix object. Each row represents a day.
#'
#' @author Azeddine Frimane \email{Azeddine.frimane@@uit.ac.ma; Azeddine.frimane@@yahoo.com}
#'
#' @examples
#'
#' # EHI generation with 1-min time step in Kenitra, Morocco.
#' extEHIcalc  <- function(phi = 34, lg = -6, tStep = 60)
#'
#' @export

#==================================================================================================
extEHIcalc <- function(phi = 0, lg = 0, tStep = 60) {
#==================================================================================================

  calc <- function(t) {

    w <- 15*(t*tStep/3600 - 12)
    delta <- 23.45*sin(0.01745329*(nj+284))
    sinh <- cos(0.01745329*delta)*cos(0.01745329*phi)*cos(0.01745329*w) +
            sin(0.01745329*phi)*sin(0.01745329*delta)
    s <- 1367*(1-0.08547009*sin(0.01745329*delta))*sinh
    if (s < 0)    s <- 0
    return(s)
  }
  H <- matrix(NA, nrow = 366, ncol = floor(24*3600/tStep))
  for (nj in 1:366) {
    S <- sapply(c(1:floor(24*3600/tStep)), function(t) calc(t))
    H[nj, ] <- S
  }
  return(H)
}
#================================================ Fin =================================================

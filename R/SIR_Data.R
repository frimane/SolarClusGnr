
#' SIRData class
#'
#' Constructor function of the objects of SIRData class needed to perform the clustering procedure of the
#' daily clearness index distrubtions. Once you create a SIRData object, you no longer need other inputs. All the
#' rest work will done automatically.
#'
#' @usage SIR_Data(Ghi_from = NULL, Ehi_from = NULL, Ehi_to = NULL, CI = NULL, Xbins = NULL)
#'
#' @param Ghi_from list of the global solar irradiance data with 10-min or 15-min loggin interval. Each component represents a day.
#' Since the day lengths throughout the year are not equal, the class list is the most suitable.
#' @param Ehi_from list of the corresponding extraterrestrial solar irradiance data of Ghi_from with
#' the same time resolution (see \code{\link{extEHIcalc}}. It can generates Ehi_from if you dont have it).
#' @param Ehi_to list of the extraterrestrial solar irradiance data of the 1-min resolution (see \code{\link{extEHIcalc}}).
#' @param CI list of the clearness index data. to be used if Ghi_from list is not available,
#'  else it is set to NULL.
#' @param Xbins numeric(1) represents the number of bins to create the matrix of the daily clearness index distributions.
#'
#' @return an object of class SIRData
#'
#' @author
#' Azeddine Frimane \email{Azeddine.frimane@@uit.ac.ma; Azeddine.frimane@@yahoo.com}
#'
#' @examples
#' data("GHI_10_min", "EHI_10_min", "EHI_1_min")
#'
#' # creat an objet of class SIRData
#' newSIRDataObj <- SIR_Data(GHI_10_min, EHI_10_min, EHI_1_min, Xbins = 13)
#'
#' @export
#============================================================================================
SIR_Data <- function(Ghi_from = NULL, Ehi_from = NULL, Ehi_to = NULL, CI = NULL, Xbins = NULL) {
#============================================================================================

  if(is.null(CI)){
  # +++ Construire la liste de Clearness Index si elle est NA +++
  CI <- mapply("/", Ghi_from, Ehi_from, SIMPLIFY = FALSE)
  CI <- lapply(CI, function(x) {
     x[ x > 1 ] <- 1
     x[ x < 0 ] <- 0
     return(x)
  })
  } else {
    # +++ Sinon construire la liste du rayonnement global  +++
    Ghi <- mapply("*", CI, Ehi_from, SIMPLIFY = FALSE)
  }

  # +++ Construire la matrice des densités +++
  dcidL <- lapply(CI, function(x) (hist(x,breaks = seq(0, 1, by = 1/Xbins),
                                          plot = FALSE)$density)/Xbins)
  dcidM<- matrix(unlist(dcidL), nrow = length(dcidL), ncol = Xbins, byrow = TRUE)

  # +++ Creer la liste de tous les donnees +++
  values <- list(clearnessIndex = CI, matrixCI = dcidM, globalRadiation = Ghi_from,
                 extraRadiation = Ehi_to)

  # +++ Creer la class +++
  attr(values, "class") <- "SIRData"
  return(values)
}
#=========================================== Fin ===============================================

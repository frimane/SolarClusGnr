
#' SIRData class
#'
#' Constructor function of the objects of SIRData class needed to perform the clustering procedure of the
#' daily clearness index distrubtions. Once you create a SIRData object, you no longer need other inputs. All the
#' rest work will done automatically.
#'
#' @usage SIR_Data(Ghi_from = NULL, Ehi_from = NULL, Ehi_to = NULL, Xbins = NULL)
#'
#' @param Ghi_from list of the global solar irradiance data with 10-min or 15-min loggin interval. Each component represents a day.
#' Since the day lengths throughout the year are not equal, the class list is the most suitable.
#' @param Ehi_from list of the corresponding extraterrestrial solar irradiance data of Ghi_from with
#' the same time resolution.
#' @param Ehi_to list of the extraterrestrial solar irradiance data of the 1-min resolution.
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
SIR_Data <- function(Ghi_from = NULL, Ehi_from = NULL, Ehi_to = NULL, Xbins = NULL)
  {
#============================================================================================

# In the name of Allah the Merciful
# =================================

  # +++ Construire la liste de Clearness Index +++
  CI <- mapply("/", Ghi_from, Ehi_from, SIMPLIFY = FALSE)
  CI <- lapply(CI, function(x) {
     x[ x > 1 ] <- 1
     x[ x < 0 ] <- 0
     return(x)
  })

  # +++ Construire la matrice des densités +++
  dcidL <- lapply(CI, function(x) (hist(x,breaks = seq(0, 1, by = 1/Xbins),
                                          plot = FALSE)$density)/Xbins)
  dcidM<- matrix(unlist(dcidL), nrow = length(dcidL), ncol = Xbins, byrow = TRUE)

  # +++ Le type du jours concernant standard deviation +++
  Tj <- lapply(CI, function(x) abs(diff(x)))
  a <- unlist(lapply(Tj, function(x) sd(x, na.rm = T)))
  b1 <- quantile(a, probs = .7, na.rm = T)
  b2 <- quantile(a, probs = .1, na.rm = T)
  b3 <- quantile(a, probs = .95, na.rm = T)

  Tjour <- mapply(function(x, j) {

    y <- numeric(length(x))
    if(j >= b1){
      if(j < b3) {
        for (i in 1:length(x)) {
          if(x[i] >= b1 || is.na(x[i])) y[i] <- 1
          else y[i] <- 0
        }
        } else {
          for (i in 1:length(x)) {
            if(x[i] >= b2 || is.na(x[i])) y[i] <- 1
            else y[i] <- 0
          }
        }
    } else {
      for (i in 1:length(x)) {
        y[i] <- 0
      }
    }
    return(c(y, 0))
  }, Tj, a, SIMPLIFY = F)

  # +++ Creer la liste de tous les donnees +++
  values <- list(clearnessIndex = CI, matrixCI = dcidM,
                 extraRadiation = Ehi_to, Tjour = Tjour)

  # +++ Creer la class +++
  attr(values, "class") <- "SIRData"
  return(values)
}
#=========================================== Fin ===============================================


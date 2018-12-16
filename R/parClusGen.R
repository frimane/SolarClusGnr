
#' genData class
#'
#' Constructor function of the objects of genData class, needed to perform the generation procedure of hight resolution
#' solar irradiance data. It is based on the resulting object dpgmmRes of the clustering procedure of the considered solar
#' irradiance data. It also can be used to generate data based on similar climate conditions by setting genNewSqClass as TRUE
#' and a dpgmmRes object corresponding to the similar climate (other than the considered solar irradiance data).
#'
#' @usage parClusGen(measures = NULL, dpgmmRes = NULL, genNewSqClass = FALSE, Xbins = NULL)
#'
#' @param measures object of class SIRData.
#' @param dpgmmRes object of class clusData.
#' @param genNewSqClass boolean variable. If TRUE the function assigne days to classes based on similar climate approach. Default is
#' FALSE.
#' @param Xbins numeric(1) used if genNewSqClass is TRUE. It represents the number of bins to create the matrix of the
#' daily clearness index distributions.
#'
#' @return an object of class genData
#'
#' @author Azeddine Frimane \email{Azeddine.frimane@@uit.ac.ma; Azeddine.frimane@@yahoo.com}
#'
#' @examples
#' data("clusData_obj")
#' data("SIRData_obj")
#'
#' # creat an object of class genData.
#' newgenDataObj <- parClusGen(SIRData_obj, clusData_obj)
#'
#' @importFrom markovchain markovchainFit
#'
#' @export
#=============================================================================================
parClusGenn <- function(measures = NULL, dpgmmRes = NULL, genNewSqClass = FALSE, Xbins = NULL)
  {
#=============================================================================================

# In the name of Allah the Merciful
# =================================

  # +++ Fonction de vraisemblance +++
  Lmultivar_t <- function (d, p.) {

    0.5*Xbins*(log(p.[[3]]/(p.[[3]]+1))-log(pi)) + lgamma(0.5*(p.[[4]]+1)) -
      lgamma(0.5*(p.[[4]]+1-Xbins)) + 0.5*p.[[4]]*log(det(p.[[2]])) -
      0.5*(p.[[4]]+1)*log(det(p.[[2]]+(p.[[3]]/(p.[[3]]+1))*tcrossprod(d-p.[[1]])))
  }

  meas <- measures$clearnessIndex

  if(genNewSqClass == TRUE){
    # +++ Assigner des classes aux donnees non labiliser en maximisant leur
    # vraisemblance +++
    dcidX <- lapply(meas, function(x) (hist(x,breaks = seq(0, 1, by = 1/Xbins), plot = FALSE)$density)/Xbins)

    lng <- dpgmmRes$cl_seq
    par <- dpgmmRes$like_par
    vraiSembl <- list()
    vraiSemblable <- matrix(nrow = length(dcidX), ncol = max(lng))
    for (i in 1:max(lng)) {
      vraiSembl[[i]] <- unlist(lapply(dcidX, function(x) Lmultivar_t(x, par[[i]])))
      vraiSemblable[ ,i] <- vraiSembl[[i]]
    }

    stat <- as.numeric(apply(vraiSemblable, 1, function(x) which(x == max(x, na.rm = T))))

  } else {

    stat <- dpgmmRes$cl_seq
  }

  # +++ COnstruire la chaine de Markov intraclasses +++

  interv <- list()
  interval <- list()
  for (j in 1:max(stat)) {
    idx <- which(stat == j)
    cilass <- list()
    for (t in 1:length(idx)) {
      cilass[[t]] <- meas[[idx[t]]]
    }
    interv[[j]] <- lapply(cilass, function(x) {
      y <- findInterval(x, seq(.01, 1, .01))
      return(as.character(y))
      }
      )

    interval[[j]] <- unlist(interv[[j]])
  }

  fitMarkov <- lapply(interval, function(x) markovchain::markovchainFit(data = x,
                                                                        sanitize = T,
                                                                        method = "laplace"))

  # +++ La liste des parametres pour la generation des donnees +++
  ret <- list(FittedMarkovChain = fitMarkov, SequenceOfStates = stat,
              Measures = measures, Tjour = measures$Tjour)

  # +++ Creer la class +++
  attr(ret, "class") <- "genData"
  return(ret)
  }
#=========================================== Fin ================================================


#' Solar irradiance data generation
#'
#' Generic function to generate high resolution solar irradiance data. It can genarate 1-min data from 10-min and
#' 15-min solar irradiance databases.
#'
#' @usage GenData(clusChar = NULL, stp_from = "10 minutes", n_cores = 1, mc = 400)
#'
#' @param clusChar object of class genData.
#' @param stp_from charachter object with values "10 minutes" or "15 minutes". It represents the original time resolution of the solar irradiance data
#' @param n_cores numreic(1), is the used number of cores. Defaut is 1. This function uses the parallel calculation for time optimisation.
#' For WINDOWS users it must be set to 1, since windows platform do not deals with fork() function.
#' @param mc numreic(1)represents the number of iterations to choose the most probable paths.
#' Optimal values around 400-600.
#'
#' @return a list object of 1-min solar irradiance data. Each element represents a day.
#'
#' @author Azeddine Frimane \email{Azeddine.frimane@@uit.ac.ma; Azeddine.frimane@@yahoo.com}
#'
#' @examples
#' data("genData_obj")
#' data("GHI_10_min")
#'
#' # data generation. For WINDOWS users n_cores must be 1.
#' GHI1min_gen <- GenData(clusChar = genData_obj, stp_from = "10 minutes", n_cores = 1)
#'
#' @importFrom markovchain markovchainSequence
#' @importFrom parallel mcmapply
#' @importFrom pbmcapply pbmcmapply
#' @importFrom utils tail
#'
#' @export
#=============================================================================================
GenData <- function(clusChar = NULL, stp_from = "10 minutes", n_cores = 1, mc = 400) {
#=============================================================================================

  ###########################################################################################
  ##### ******************************* Initialisation  ******************************* #####
  ###########################################################################################

  # +++ Determiner le nombre de fois du lancement des boucles d'echantillonage +++
  # +++ Dans ces boucle le nombre d'echantillons est stp +++
  if(stp_from == "15 minutes"){

      stp <- 15
  } else if(stp_from == "10 minutes"){

      stp <- 10
  } else {

    print("Insert stp_from?")
  }

  Meas <- clusChar$Measures
  # +++ Les intervales ou les pas d'echantillonnage de CI +++
  intervals <- lapply(Meas$clearnessIndex, function(x) {
    y <- findInterval(x, seq(.01, 1, .01))
    return(as.character(y))
  }
  )

  # +++ La sequence des classes (differencie par S majuscule
  # par rapport au states noms des colonnes des matrices de transition) +++
  State <- clusChar$SequenceOfStates

  # +++ Les matrice de transition pour chaque cluster +++
  Trans <- clusChar$FittedMarkovChain

  # +++ Detriminer les tiks pour l'echantillonage de CI +++
  Tiks <- seq(.01, 1, .01) - .005

  ###########################################################################################
  ##### ********************************* Main function ******************************* #####
  ###########################################################################################

  Generation <- function(x, i, stp){

    # # +++ La matrice qui contient les intervale echantilloner pour chaque jour +++
    # # +++ le probleme est modeliser de tel sorte que chaque ligne de cette matrice
    # # contient les intervales echantilloner pour chaque 10 echantillons +++
    c.name <- matrix(nrow = length(x), ncol = stp)

    c.name[1,1] <- x[1]

    # +++ Markov chain Object, la sequence des etats et la matrice de transition +++
    te <- Trans[[State[i]]]$estimate
    st <- attributes(te)$states
    tm <- attributes(te)$transitionMatrix

    for(j in 1:length(x)) {

      # +++ Verifier si la matrice des transition contient l'etat initial, sinon
      # mettre a sa place l'etat la plus prbable  +++
      if(is.na(x[j])){
        c.name[j, ] <- rep(NA, stp)
      } else {

        if(x[j] %in% st == FALSE) {
        for (l in j:length(x)) {
          repeat{
            x[j] <- x[l+1]
            if(l+1 == length(x)){
              x[j] <- sample(st, 1)
            } else if(x[j] %in% st == TRUE) {
              break
            } else {
              x[j] <- sample(st, 1)
            }
          }
        }
      }

      # +++ Assurer la continuité des etats +++
      xj1 <- x[j+1]
      if(is.na(xj1)){
        xj1 <- sample(st, 1)
      } else if(xj1 %in% st == FALSE) {
        for (l in j:length(x)) {
          repeat{
            xj1 = x[l+1]
            if(xj1 %in% st == TRUE) { break
            } else {
              xj1 = sample(st, 1)
            }
          }
        }
      }

      # +++ Monter-carlo pour choisir les etats +++
      re.c.name <- replicate(mc, markovchain::markovchainSequence(n = stp, markovchain = te,
                                                                  t0 = x[j], include.t0 = F))

      # +++ Choisir les colonnes qui finissent par l'etat qui vient +++
      fin.name <- as.matrix(re.c.name[, tail(re.c.name, 1) == xj1])

      # +++ Choisir le chemin le plus probable +++
      if(ncol(fin.name) != 0) {

        pb.c.name <- as.numeric(apply(fin.name, 2, function(nm) {
          y <- rep(1, length(nm))
          for (k in 2:length(nm)) {

            y[k] <- -log(tm[nm[k-1], nm[k]])
          }
          a <- sum(y)
          ab <- sd(as.numeric(nm))
          abc <- a*ab
          return(abc)
        }
        ))

        c.name[j, ] <- fin.name[ , which(pb.c.name == min(pb.c.name))[1]]

        # +++ Sinon repliquer NA values +++
      } else {
        c.name[j, ] <- rep(NA, stp)
      }
      }
    }

    return(c(t(c.name)))
  }

  ###########################################################################################
  ##### ******************************** Data Generation ****************************** #####
  ###########################################################################################

  # +++ Generer les donnees +++
    G.data <- pbmcapply::pbmcmapply(Generation, intervals, seq_along(intervals), stp, SIMPLIFY = F,
                       mc.cores = n_cores)

    # +++ la liste qui contient les echantillons +++
    rCI <- lapply(G.data, function(x){
      y <- numeric(length(x))
      for (i in 1:length(x)) {
        y[i] <- Tiks[as.numeric(x[i]) + 1] + runif(1,min = -.003, max = .003)
      }
      return(y)
    })

 GHi <- mapply("*", rCI, Meas$extraRadiation,SIMPLIFY = F)

return(GHi)
}
#=========================================== Fin ================================================



#' DPGMMclus S3 Method
#'
#' Nonparametric Bayesian Dirichlet-Gaussian clustering of daily clearness index distributions. It can be also used to perform any
#' data clustering of class matrix other than irradiance data of class SIRData.
#'
#' @usage DPGMMclus(obj, n.iter, n.burn, ci, alfa)
#'
#' @param obj object of class SIRData (see \code{\link{SIR_Data}}) or a matrix object containing the clearness index distributions.
#' @param n.iter numeric(1) represents the number of iterations
#' @param n.burn numeric(1) represents the number of burn-in iterations (ignored iterations)
#' @param ci numeric vector that contain the initialization of the indicator variables (the initial assignement of data to clusters).
#' @param alfa numeric(1) represents the concentration parameter. Default is 0.5.
#'
#' @return an object of clusData class containing:
#' \item{cl_seq}{numeric vector represents the class sequence.}
#' \item{like_par}{list0object represents the parameters of each class. 1st element is the mean (numeric vector), 5th element is the precision matrix
#' (inverse of the co-variance matrix) of the class.}
#' \item{lik_con_par}{numeric(1), represents the inferred concentration parameter alf.}
#' \item{Con_par_sam_seq}{numeric vector, represents the posterior distribution of alf.}
#' \item{seq_cl_num}{numeric vector, represents the distribution of the class numbers.}
#' \item{bet_seq}{numeric vector, represents the destribution of the beta parametrer of the precision matrix.}
#' \item{like_cl}{numeric vector, represents the number of elements of each class.}
#' \item{calc_time}{numeric(1) represents the computing time consumed.}
#'
#' @author Azeddine Frimane \email{Azeddine.frimane@@uit.ac.ma; Azeddine.frimane@@yahoo.com}
#'
#' @examples
#' data("SIRData_obj")
#'
#' newClustering <- DPGMMclus(SIRData_obj, n.iter = 1000, n.burn = 500, ci = 5, alfa = 1)
#' # for class ploting see \code{\link{clPlot}}
#'
#' @export
#============================================================================================
DPGMMclus <- function(obj, n.iter, n.burn, ci, alfa) {

  UseMethod("DPGMMclus", obj)
}
#============================================================================================
#' @export
#============================================================================================
DPGMMclus.SIRData <- function(obj, n.iter = 1000, n.burn = 500,
                             ci = 30, alfa = .5) {

  ClusDirGauss(meas = obj$matrixCI, n_iter = n.iter, n_burn = n.burn,
               zi = sample(1:ci, nrow(obj$matrixCI), replace = T), alf = alfa)
}
#============================================================================================
#' @export
#============================================================================================
DPGMMclus.default <- function(obj, n.iter = 2000, n.burn = 500,
                              ci = 30, alfa = .5) {

  ClusDirGauss(meas = obj, n_iter = n.iter, n_burn = n.burn,
               zi = sample(1:ci, nrow(obj), replace = T), alf = alfa)
}
#=========================================== Fin ================================================

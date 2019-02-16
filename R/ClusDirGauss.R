
#' Dirichlet process Gaussian Mixture Model -- DPGMM --
#'
#' Generic function of the nonparametric Bayesian Dirichlet-Gaussian clustering.
#'
#' @importFrom Matrix nearPD
#' @importFrom matrixcalc is.positive.definite
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats pnorm
#' @importFrom stats rWishart
#' @importFrom stats rgamma
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom stats var
#'
#' @author Azeddine Frimane \email{Azeddine.frimane@@uit.ac.ma; Azeddine.frimane@@yahoo.com}
#'
#' @export
#'
#============================================================================================
ClusDirGauss <- function(meas , n_iter = NULL, n_burn = NULL, zi = NULL, alf = NULL) {
#============================================================================================

#************************
# all needed functions  *
#************************

# Update component parameters
updateParameters <- function (d, n) {

  if (is.vector(d)) d <- matrix(d, nrow = 1)
  r <- n + h[[3]]
  b <- n + h[[4]]
  u <- (h[[3]]*h[[1]] + colSums(d))/r
  t <- h[[2]] - r*tcrossprod(u) + h[[3]]*tcrossprod(h[[1]]) +
    Reduce('+', lapply(c(1:nrow(d)), function(x) tcrossprod(d[x, ])))
  if(is.positive.definite(t) == FALSE) t <- as.matrix(nearPD(t)$mat)
  s <- matrix(rWishart(1, b, chol2inv(chol(t))), D, D)
  m <- as.numeric(mvtnorm::rmvnorm(1, mean = u, sigma = chol2inv(chol(r*s))))
  return(list(m, t, r, b, s))
}

# Update hyperparameters
updateHyperparameters <- function() {

  y <- Reduce('+',lapply(p, function(x) crossprod(x[[1]]-h[[1]],x[[5]])%*%(x[[1]]-h[[1]])))
  y. <- Reduce('+',lapply(p, function(x) x[[5]]))
  y.. <- h[[3]]*Reduce('+',lapply(p, function(x) x[[5]]%*%x[[1]]))
  s <- solve(h[[3]]*y. + dpr)
  b <- McMcb(updateBeta, h[[4]], 1)
  m <- as.numeric(rmvnorm(1, mean=as.numeric(s%*%(y..+dpr%*%dav)), sigma=s))
  t <- b*matrix(rWishart(1, D+max(zi)*h[[4]], solve(h[[4]]*y.+dpr)), D, D)
  r <- rgamma(1, shape = 0.5*(1+max(zi)), rate = 0.5*(1+y))

  return(list(m, t, r, b))
}

# update the DP concentration parameter alfa.
updateAlfa <- function(x) {

  return((max(zi)-0.5)*x - 0.5*exp(-x) + lgamma(exp(x)) - lgamma(exp(x)+ N))
}

# update the beta hyperparameters.
updateBeta <- function(x){
  x - max(zi)*Reduce('+' ,lapply(c(1:D), function(y) lgamma(0.5*(exp(x)+y-D)))) +
    0.5*D*max(zi)*exp(x)*(x-log(2)) - (0.5*D)/(exp(x)-D+1) - 1.5*log(exp(x)-D+1) +
    0.5*exp(x)*Reduce('+',lapply(p, function(y) log(det(y[[5]])*det(h[[2]]/h[[4]])) -
                                   sum(diag(y[[5]]%*%h[[2]]/h[[4]]))))
}

# markov chain monte carlo for sampling alfa (log-normal function)
McMc <- function(target, x, sd) {

  target.x.current <- target(log(x))
  x.current <- log(x)
  x.proposed <- rnorm(1, x.current, sd)
  target.x.proposed <- target(x.proposed)
  log.acceptance <- target.x.proposed - target.x.current

  r <- runif(1)
  if (r < exp(log.acceptance)) {
    x.current <- x.proposed
  }
  return(exp(x.current))
}

McMcb <- function(target, x, sd ) {

  target.x.current <- target(log(x))
  x.current <- log(x)

  repeat {
    x.proposed <- rnorm(1, x.current, sd)
    if(x.proposed > log(D-1)) break
  }

  target.x.proposed <- target(x.proposed)
  log.acceptance <- target.x.proposed - target.x.current
  log.acceptance <- log.acceptance + pnorm(x.current, log.p = T) - pnorm(x.proposed, log.p = T)

  r <- runif(1)
  if (r < exp(log.acceptance)) {
    x.current <- x.proposed
  }
  return(exp(x.current))
}

# posterior predictive ~ multivariate t-distribution
Lmultivar_t <- function (d, p.) {

  0.5*D*(log(p.[[3]]/(p.[[3]]+1))-log(pi)) + lgamma(0.5*(p.[[4]]+1)) -
    lgamma(0.5*(p.[[4]]+1-D)) + 0.5*p.[[4]]*log(det(p.[[2]])) -
    0.5*(p.[[4]]+1)*log(det(p.[[2]]+(p.[[3]]/(p.[[3]]+1))*tcrossprod(d-p.[[1]])))
}

# multinomial-distribution
updateIndicator <- function () {

  pb <- c(log(n_el) + sapply(p, function(x) Lmultivar_t(meas[j, ], x)),
          log(alf) + Lmultivar_t(meas[j, ], h))
  pb <- exp(pb - max(pb))
  pb <- pb/sum(pb)

  x <- runif(1)
  x. <- 0
  for (k in 1:length(pb) ) {
    x. <- x. + pb[k]
    if (x < x.) break;
  }
  return(k)
}

#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

# MCMC parameters
# ---------------

likely_max <- -Inf
n_cl <- c()
d_alf <- c()
d_beta <- c()
likely_cluster <- 0

# data
# ----

N <- nrow(meas)
D <- ncol(meas)

# hyperpriors
# -----------
dav <- apply(meas,2,mean)
dco <- as.matrix(nearPD(var(meas))$mat) #--------------------------- not.positive.definite
dpr <- solve(dco)#chol2inv(chol(dco))

# initial cluster structure
# -------------------------
n_el <- as.numeric(table(factor(zi)))
B <- D
h <- list(as.vector(mvtnorm::rmvnorm(1, dav, dco)),
          B*matrix(rWishart(1, D, dco/D), D, D),
          rgamma(1, 1, 1),
          B)

set.seed(19)
#--------------------------------------------------------------------------------------------------#
t1 <- Sys.time()
#--------------------------------------------------------------------------------------------------#

for (i in 1:n_iter) {

  p <- lapply(c(1:max(zi)), function(x) updateParameters(meas[zi==x,], n_el[x]))
  h <- updateHyperparameters()

  for (j in 1:N) {

    n_el[zi[j]] <- n_el[zi[j]] - 1

    if (n_el[zi[j]] == 0) { #--------------------------- relabel clusters
      p[[zi[j]]] <- NULL
      n_el <- n_el[-zi[j]]
      zi[zi>zi[j]] <- zi[zi>zi[j]]-1
    }

    zi[j] <- updateIndicator() # ------------------------------- drawing a new indicator variable

    if (length(n_el) < zi[j]) {
      n_el[zi[j]] <- 1
      p[[zi[j]]] <- updateParameters(rep(0,D), 1)
    } else { n_el[zi[j]] <- n_el[zi[j]] + 1 }
  }

  if (i > n_burn) {

    likely <- sum(sapply(c(1:N), function(x)

                           (h[[4]]-D + 1)*log(det(p[[zi[x]]][[5]])) - h[[4]]*sum(diag(p[[zi[x]]][[5]]%*%h[[2]]/h[[4]])) -
                           c(h[[3]]*crossprod((p[[zi[x]]][[1]]-h[[1]]), p[[zi[x]]][[5]]%*%(p[[zi[x]]][[1]]-h[[1]]))) -
                           c(crossprod((meas[x,]-p[[zi[x]]][[1]]),p[[zi[x]]][[5]]%*%(meas[x,]-p[[zi[x]]][[1]]))) ))

    if (likely > likely_max) {

      likely_max <- likely
      likely_cluster <- zi
      likely_par <- p
      likely_iter <- i
      likely_alf <- alf
    }

  }

  d_alf <- c(d_alf, alf)
  d_beta <- c(d_beta, h[[4]])
  n_cl <- c(n_cl, max(zi))
  alf <- McMc(updateAlfa, alf, 1)

  cat("\nIteration:",i,"","", "number_of_classes:",max(zi),":", table(zi),"","", "likely_number_of_class:",max(likely_cluster),"\n")
}

numberLatentClasses <- table(factor(likely_cluster))
#****************************************************
t2 <- Sys.time()
howMuchTime <- difftime(t2, t1)

res <- list(cl_seq= likely_cluster, like_par=likely_par, lik_con_par=likely_alf, Con_par_sam_seq=d_alf,
            seq_cl_num=n_cl, bet_seq=d_beta,
            like_cl=numberLatentClasses, calc_time=howMuchTime)

# +++ Creer la class +++
attr(res, "class") <- "clusData"
return(res)
}
#=========================================== Fin ===============================================


#' class plot
#'
#' Generic function to generate plot of classes.
#'
#' @usage clPlot(meas = NULL, dpgmm = NULL, nc = 2)
#'
#' @param meas object of class SIRData.
#' @param dpgmm object of class clusData.
#' @param nc number of column of plot.
#'
#' @return Multiplot panel. Each plot represents a class. The histograms represent the mean
#' clearness Index distribution of the class and the bars represent the standard deviation
#' of the class.
#'
#' @author Azeddine Frimane \email{Azeddine.frimane@@uit.ac.ma; Azeddine.frimane@@yahoo.com}
#'
#' @examples
#'
#' data("SIRData_obj")
#' data("clusData_obj")
#'
#' # three classes of GHI.
#' clPlot(SIRData_obj, clusData_obj, 2)
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_errorbar
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 margin
#' @importFrom graphics hist
#'
#' @export
#============================================================================================
clPlot <- function(meas = NULL, dpgmm = NULL, nc = 2) {
#============================================================================================

  dcid <- meas$matrixCI
  b <- ncol(dcid)
  n <- length(dpgmm$like_par)
  mn <- sapply(1:n, function(i) apply(dcid[dpgmm$cl_seq == i,], 2, mean))
  SD <- sapply(1:n, function(i) apply(dcid[dpgmm$cl_seq == i,], 2, sd))
  df <- as.data.frame(cbind(x = seq(1/(2*b),(2*b-1)/(2*b),1/b), mn, SD))

  pl <- list()
  for (i in 1:(.5*(ncol(df)-1))) {
      p <- eval(substitute(
                ggplot(df, aes(x = df[,1], y = df[,i+1])) +
                geom_errorbar(aes(ymin = df[,i+1], ymax = df[,i+1] + df[,i+(.5*(ncol(df)-1))+1]),
                              colour="grey40", width=.03, size=0.2) +
                geom_bar(stat="identity",fill=I("blue"), col=I("white")) +
                ylab("Density") + xlab("" ~ k[t]) +
                theme(panel.background = element_rect(fill = 'gray97'), title = element_text(size=4.5),
                     axis.text.y = element_text(angle = 90, hjust = 1),
                     axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
                    axis.title = element_text(size=10.5))
      ,list(i = i)))
      pl[[i]] <- p
    }

  multiplot(plotlist = pl, cols = nc)
}
#===================================== Fin =================================================

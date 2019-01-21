#' The Generic \code{plot} Function for Object of \code{BootstrapCI} Class
#'
#' @name plot.BootstrapCI
#' @aliases plot.bootstrapCI
#' @param cutpoint an object of \code{BootstrapCI} class, returned from Function \code{BootstrapCI}.
#' @param main The title of the plot.
#' @param cex.text The font size for labeling the four sets of bootstrap CIs.
#' @param a0,b0,gap Figure parameters for setting some margins in the plot to make it look nicers.
#' No need to adjust in general.
#' @details
#' The (generic) plot method for an \code{BootstrapCI} object. It plots the histogram of the bootstrapped statistics,
#' together with four sets of bootstrap confidence intervals.
#'
#' @return Plots of the results from \code{BootstrapCI}.
#' @examples
#'       data(fertility)
#'       form0 <- diagnosis ~  factor(season) + age + accident + fever + alcohol + sitting
#'       out.boots <- BootstrapLogitBestcut(formula=form0, dat=fertility, B=100,
#'                  trace=FALSE, conf.level=0.95, CV.method="none",  method.OptimalCutpoints = "Youden")
#'       out <- BootstrapCI(out.boots$c.bootstrap, out.boots$N.Bn, t.0=out.boots$c.0, trace=FALSE,
#'                  conf.level=0.95, print.it=TRUE)
#'       out
#'       plot(out)
#'
#' @seealso \code{\link{BootstrapLogitBestcut}}, \code{\link{BootstrapCI}}
#' @references
#'\itemize{
#' \item Zhang, Z., Su, X., et al. (2019+). Bootstrap Confidence Intervals for Optimal Cutoff Point
#' in Logistic Regression. Submitted.
#' }
#' @import OptimalCutpoints
#' @export


plot.BootstrapCI <- function(cutpoint, main="Histogram of Bootstrapped Optimal Cutpoints",
                    cex.text=0.8, a0=5, b0=a0-1+.2, gap=0.2){
  c.bootstrap <- cutpoint$t.bootstrap
  hist.out <- hist(c.bootstrap, nclass=50, plot =FALSE)
  # hist.out$counts <- hist.out$counts/sum(hist.out$counts)
  max.p <- max(hist.out$counts); min.p <- min(hist.out$counts)
  range.p <- max.p - min.p
  plot(hist.out, ylim=c(-range.p/a0, max.p), xlim=c(0,1), col="orange",
       main=main, xlab="", ylab="Frequency", xaxt='n', yaxt='n')
  mtext(text="Bootstrapped Cutpoints", side=1, line=2.5)
  # abline(h=0)
  # text(x=0:10/10, y=-range.p/20, labels =0:10/10)
  axis(2, at=seq(0, max.p, length.out=5), labels=ceiling(seq(0, max.p, length.out=5)))
  axis(1, at=0:10/10, labels=0:10/10)
  # lines(density(c.bootstrap, adjust=2), lty="dotted", col="orange")
  rug(jitter(c.bootstrap), ticksize = 0.01, side=1, col="gray50")
  # abline(v=cutpoint$result["t.0"], lwd=2, col="deepskyblue4")
  # abline(v=cutpoint$result["t.bagged"], lwd=2, col="darkolivegreen4")

  # PERCETILE CI
  x.percentile <- (-range.p/a0)/b0
  LB.percentile <- as.numeric(cutpoint$result["LB2.bootstrap.2.5%"])
  UB.percentile <- as.numeric(cutpoint$result["UB2.bootstrap.97.5%"])
  arrows(x0=LB.percentile, y0= x.percentile, x1=UB.percentile, y1=x.percentile,
         length=0.05, angle=90, lwd=2, code=3, col="blue") # PERCENTILE CI
  text(x=UB.percentile+gap-0.05, y=x.percentile, labels="Percentile", cex=cex.text, col="blue")

  # BOOTSTRAP CI
  x.bootstrap <- (-range.p/a0)*2/b0
  LB.bootstrap <- cutpoint$result["LB1.bootstrap"]; UB.bootstrap <- cutpoint$result["UB1.bootstrap"]
  arrows(x0=LB.bootstrap, y0=x.bootstrap, x1=UB.bootstrap, y1=x.bootstrap,
         length=0.05, angle=90, lwd=2, code=3, col="blue") # BOOTSTRAP CI
  points(x=cutpoint$result["t.0"], y=x.bootstrap, pch=19, cex=1.2, col="blue")
  text(x=UB.bootstrap+gap-0.05, y=x.bootstrap, labels="Bootstrap", cex=cex.text, col="blue")

  # IJ-BASED CI
  x.IJ <- (-range.p/a0)*3/b0
  LB.IJ <- cutpoint$result["LB1.IJ"]; UB.IJ <- cutpoint$result["UB1.IJ"]
  arrows(x0=LB.IJ, y0=x.IJ, x1=UB.IJ, y1=x.IJ, length=0.05, angle=90, lwd=2, code=3, col="green4") # BOOTSTRAP CI
  points(x=cutpoint$result["t.bagged"], y=x.IJ, pch=19, cex=1.2, col="green4")
  text(x=UB.IJ+gap-0.1, y=x.IJ, labels="IJ", cex=cex.text, col="green4")

  # BIAS-CORRECTED IJ-BASED CI
  x.IJ.unbiased <- (-range.p/a0)*4/b0
  LB.IJ.unbiased <- cutpoint$result["LB2.IJ"]; UB.IJ.unbiased <- cutpoint$result["UB2.IJ"]
  arrows(x0=LB.IJ.unbiased, y0=x.IJ.unbiased, x1=UB.IJ.unbiased, y1=x.IJ.unbiased, length=0.05, angle=90, lwd=2, code=3, col="green4") # BOOTSTRAP CI
  points(x=cutpoint$result["t.bagged"], y=x.IJ.unbiased, pch=19, cex=1.2, col="green4")
  text(x=UB.IJ.unbiased++gap, y=x.IJ.unbiased, labels="Bias-Corrected IJ", cex=cex.text, col="green4")
}



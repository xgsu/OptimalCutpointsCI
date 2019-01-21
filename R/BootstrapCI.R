#' Extract Four Sets of Confidence Intervals (CI) from the Bootstrap Results
#'
#' @name BootstrapCI
#' @aliases bootstrapCI, bootsCI, BootsCI
#' @param t.bootstrap A vector of bootstrapped statistics of dimention \eqn{B}.
#' @param N.Bn  The incidence matrix of dim \eqn{latex}{B x n} associated with the bootstrap
#' resampling procedure, where element \eqn{latex}{B_{bi}} is the frequency that the \eqn{i}-th observation
#' shows up in the \eqn{b}-th bootstrap sample.
#' @param t.0 The statistics computed on basis of the original sample.
#' @param trace A logical value. If \code{TRUE}, information on progress is shown. The default is \code{FALSE}.
#' @param conf.level A numerical value with the confidence level for the construction of the
#' confidence intervals. The default value is 0.95.
#' @param print.it A logical value. If \code{TRUE}, a summary of results is printed. Defaulted as \code{TRUE}.
#' @details
#'       Four sets of confidence intervals are extracted: the bootstrap confidence interval and the percentile CI for
#'       the sample estimate of the optimal cutpoint; the infinitesimal jackknife (IJ) confidence interval and its
#'       bias-corrected version for the bagging estimator of the optimal cutpoint. The bagging estimator is more
#'       precise with smaller variation. The bias correction is helpful when \eqn{B} is small.
#'
#'       This function is used inside of Function \code{LogitOptimalCutpointsCI} in this package, but it is a generaic function that
#'       can be used for other general purposes.
#' @return
#' The following results are returned:
#' \describe{
#' \item{result}{A vector of results including the sample estimate \code{t.0}, its bootstrap SE \code{SE.bootstrap},
#'       the lower/upper bound of the bootstrap CI (\code{LB1.bootstrap}, \code{UB1.bootstrap}), the bagging
#'       estimate \code{t.bagged}, the IJ-based SE \code{SE.IJ}, the lower/upper bound of the IJ-based CI
#'       (\code{LB1.IJ}, \code{UB1.IJ}), the bias estimate \code{bias} for variance, the bias-corrected IJ variance
#'       estimate \code{Vc.IJ}, the bias-corrected IJ SE \code{SEc.IJ}, the lower/upper bound of the bias-corrected
#'       IJ confidence interval (\code{LB2.IJ}, \code{UB2.IJ}).}
#' \item{t.bootstrap}{A copy of the \eqn{B}-dimensional vector of bootstrapped statistics;}
#' \item{N.Bn}{A copy of the bootstrap incidence matrix;}
#' \item{V.IJ}{The IJ estimate of the variance;}
#' \item{bias}{The (upward) bias estimate for \code{V.IJ};}
#' \item{Vc.IJ}{The bias-corrected IJ variance estimate;}
#' \item{conf.level}{Confidence level used;}
#' \item{B}{The number of bootstrap samples used.}
#' }
#'
#' @examples
#'       data(fertility)
#'       form0 <- diagnosis ~  factor(season) + age + accident + fever + alcohol + sitting
#'       out.boots <- BootstrapLogitBestcut(formula=form0, dat=fertility, B=100,
#'                  trace=FALSE, conf.level=0.95, CV.method="none",  method.OptimalCutpoints = "Youden")
#'       out <- BootstrapCI(out.boots$c.bootstrap, out.boots$N.Bn, t.0=out.boots$c.0, trace=FALSE,
#'                  conf.level=0.95, print.it=TRUE)
#'       out
#'       plot.BootstrapCI(out)
#'
#' @seealso \code{\link{BootstrapLogitBestcut}}, \code{\link{OptimalThreshold}}, \code{\link{plot.BootstrapCI}}
#' @references
#'\itemize{
#' \item Zhang, Z., Su, X., et al. (2019+). Bootstrap Confidence Intervals for Optimal Cutoff Point
#' in Logistic Regression. Submitted.
#' }
#' @import OptimalCutpoints
#' @export

BootstrapCI <- function(t.bootstrap, N.Bn, t.0=NULL, trace=FALSE, conf.level=0.95, print.it=TRUE)
{
  B <- NROW(N.Bn); n <- NCOL(N.Bn)
  t.bagged <- mean(t.bootstrap, na.rm =TRUE)
  V.IJ <- sum((apply(N.Bn, 2, FUN=cov, y=t.bootstrap, method="pearson"))^2)
  bias <- var(t.bootstrap, na.rm =TRUE)*(B-1)*(n-1)/(B^2)
  Vc.IJ <- V.IJ - bias; if (Vc.IJ <0) print("Oh, No! You got negative variance after Bais Correction!")
  SE.IJ <- sqrt(V.IJ)
  SEc.IJ <- ifelse(Vc.IJ >0, sqrt(Vc.IJ), NA);

  # BOOTSTRAP CONFIDENCE INTERVALS
  # #################################
  z0 <- qnorm(p=1-(1-conf.level)/2)
  # BOOTSTRAP CIs FOR BEST CUTOFF: METHOD I
  # ------------------------------------------
  if (is.null(t.0)) t.0 <- t.bagged
  SE.bootstrap <- sd(t.bootstrap, na.rm =TRUE)
  UB1.bootstrap <- t.0 + z0*SE.bootstrap
  LB1.bootstrap <- t.0 - z0*SE.bootstrap
  CI1.bootstrap <- c(LB1.bootstrap, UB1.bootstrap)
  # BOOTSTRAP CIs FOR BEST CUTOFF: METHOD II
  # -----------------------------------------
  CI2 <- quantile(t.bootstrap, probs =c((1-conf.level)/2, 1-(1-conf.level)/2));
  LB2.bootstrap <- CI2[1]; UB2.bootstrap <- CI2[2];

  # IJ CONFIDENCE INTERVALS
  # ########################
  # WITHOUT BIAS CORRECTION
  # ------------------------
  LB1.IJ <- t.bagged - z0*SE.IJ; UB1.IJ <- t.bagged + z0*SE.IJ
  CI1.IJ <- c(LB1.IJ, UB1.IJ)
  # BIAS-CORRECTED CI
  # -------------------
  LB2.IJ <- t.bagged - z0*SEc.IJ; UB2.IJ <- t.bagged + z0*SEc.IJ
  CI2.IJ <- c(LB2.IJ, UB2.IJ)

  # PRINT OUT THE RESULTS
  if (print.it || trace) {
    cat("------------------------------------------------------\n")
    cat("The best cutoff is ", t.0, ";\n", sep="")
    cat("\t with ", conf.level, " Bootstrap CI: (", LB1.bootstrap, ", ", UB1.bootstrap, ") \n", sep="")
    cat("\t and ", conf.level, " Bootstrap Quantile CI: (",  LB2.bootstrap, ", ", UB2.bootstrap, ").\n", sep="")
    cat("------------------------------------------------------\n")
    cat("The bagged cutoff estimate is ", t.bagged, "\n")
    cat("\t with ", conf.level, " Uncorrected IJ CI: (", LB1.IJ, ", ", UB1.IJ, ");\n", sep="")
    cat("\t and ", conf.level, " Bias-Corrected IJ: (", LB2.IJ, ", ", UB2.IJ, ").\n", sep="")
    cat("------------------------------------------------------\n")
    cat("A total of B=", B, " bootstrap samples are taken.\n\n", sep="")
  }
  # OUTPUT THE RESUTLS
  result <- c(t.0=t.0, SE.bootstrap=SE.bootstrap, LB1.bootstrap=LB1.bootstrap, UB1.bootstrap=UB1.bootstrap,
              LB2.bootstrap=LB2.bootstrap, UB2.bootstrap=UB2.bootstrap,
              t.bagged=t.bagged, SE.IJ=SE.IJ, LB1.IJ=LB1.IJ, UB1.IJ=UB1.IJ,
              bias=bias, Vc.IJ=Vc.IJ, SEc.IJ=SEc.IJ, LB2.IJ=LB2.IJ, UB2.IJ=UB2.IJ)
  out <- list(result=result, t.bootstrap=t.bootstrap, N.Bn=N.Bn, V.IJ=V.IJ, bias=bias, Vc.IJ=Vc.IJ,
              conf.level=conf.level, B=B)
  class(out) <- "BootstrapCI"
  return(out)
}

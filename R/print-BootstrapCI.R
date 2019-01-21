#' The Generic \code{print} Function for Object of \code{BootstrapCI} Class
#'
#' @name print.BootstrapCI
#' @aliases print.bootstrapCI
#' @param cutpoint an object of \code{BootstrapCI} class, returned from Function \code{BootstrapCI}.
#' @examples
#'       data(fertility)
#'       form0 <- diagnosis ~  factor(season) + age + accident + fever + alcohol + sitting
#'       out.boots <- BootstrapLogitBestcut(formula=form0, dat=fertility, B=100,
#'                  trace=FALSE, conf.level=0.95, CV.method="none",  method.OptimalCutpoints = "Youden")
#'       out <- BootstrapCI(out.boots$c.bootstrap, out.boots$N.Bn, t.0=out.boots$c.0, trace=FALSE,
#'                  conf.level=0.95, print.it=TRUE)
#'       print(out)
#'
#' @seealso \code{\link{BootstrapLogitBestcut}}, \code{\link{BootstrapCI}}
#' @references
#'\itemize{
#' \item Zhang, Z., Su, X., et al. (2019+). Bootstrap Confidence Intervals for Optimal Cutoff Point
#' in Logistic Regression. Submitted.
#' }
#' @import OptimalCutpoints
#' @export


print.BootstrapCI <- function(cutpoint)
{
  result <- cutpoint$result
  t.0 <- result["t.0"]; LB1.bootstrap <- result["LB1.bootstrap"]; UB1.bootstrap <- result["UB1.bootstrap"];
  LB2.bootstrap <- result["LB2.bootstrap"]; UB2.bootstrap <- result["UB2.bootstrap"];
  LB1.IJ <- result["LB1.IJ"]; UB1.IJ <- result["UB1.IJ"]
  LB2.IJ <- result["LB2.IJ"]; UB2.IJ <- result["UB2.IJ"]
  cat("------------------------------------------------------\n")
  cat("The best cutoff is ", t.0, ";\n", sep="")
  cat("\t with ", cutpoint$conf.level, " Bootstrap CI: (", LB1.bootstrap, ", ", UB1.bootstrap, ") \n", sep="")
  cat("\t and ", cutpoint$conf.level, " Bootstrap Quantile CI: (",  LB2.bootstrap, ", ", UB2.bootstrap, ").\n", sep="")
  cat("------------------------------------------------------\n")
  cat("The bagged cutoff estimate is ", cutpoint$t.bagged, "\n")
  cat("\t with ", cutpoint$conf.level, " Uncorrected IJ CI: (", LB1.IJ, ", ", UB1.IJ, ");\n", sep="")
  cat("\t and ", cutpoint$conf.level, " Bias-Corrected IJ: (", LB2.IJ, ", ", UB2.IJ, ").\n", sep="")
  cat("------------------------------------------------------\n")
  cat("A total of B=", cutpoint$B, " bootstrap samples are taken.\n\n", sep="")
}

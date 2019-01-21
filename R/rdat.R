#' Generate some data from a logistic regression model
#'
#' @name rdat
#' @aliases rdata
#' @param n The sample size.
#' @param beta The vector of regression coefficients. The first component is intercept.
#' The number of predictors is automatically determined by the dimension of beta
#' @param x.dist Predictors are independently generated from either \code{"uniform"} or \code{"normal"}
#' distributions. Correlated predictors could have been easily used, but we didn't bother to do so.
#' @details
#' The function simualtes data from a true logistic regression model. The empirical prevalence rate for
#' \code{y=1} is printed out.
#'
#' @return A data frame \code{dat} where the first column is the binary response \code{y}, followed by
#' a number of predictors \code{x1}, \code{x2}, \eqn{\ldots}.
#' @examples
#'       beta0 <- c(-2, 2, -2, 2, -2, -2)
#'       dat <- rdat(n=1000, beta=beta0)
#'       dim(dat); head(dat)
#'
#' @seealso \code{\link{rnorm}}, \code{\link{runif}}
#' @references
#'\itemize{
#' \item Zhang, Z., Su, X., et al. (2019+). Bootstrap Confidence Intervals for Optimal Cutoff Point
#' in Logistic Regression. Submitted.
#' }
#' @import stats
#' @export


rdat <- function(n, beta, x.dist="uniform"){
	p <- length(beta)-1
	if (x.dist=="normal") x <- rnorm(n*p, 0, 1)  # MAY ADD CORRELATION
	else if (x.dist=="uniform") x <- runif(n*p, 0, 1)  # MAY ADD CORRELATION
	else stop("Either normal or uniform distributon for predictors plz.")
	X <- matrix(x, nrow=n, ncol=p)
	eta <- X%*%(beta[-1]) + beta[1]
	p0 <- as.vector(exp(eta)/(1+exp(eta)))
	y <- rbinom(n, size=1, prob=p0);
	print(mean(y))
 	dat <- data.frame(cbind(y, X))
	colnames(dat) <- c("y", paste("x", 1:p, sep=""))
	dat
}

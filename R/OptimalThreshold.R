#' A function modified from Function \code{optimal.cutpoints} in the R Package \bold{OptimalCutpoints}
#'
#' @name OptimalThreshold
#' @aliases optimal.threshold
#' @param X Either a character string with the name of the diagnostic test variable (then method
#' 'optimal.cutpoints.default' is called), or a formula (then method 'optimal.cutpoints.formula'
#' is called). When 'X' is a formula, it must be an object of class "formula". Right side of \code{~}
#' must contain the name of the variable that distinguishes healthy from diseased individuals, and
#' left side of \code{~} must contain the name of the diagnostic test variable.
#' @param status A character string with the name of the variable that distinguishes healthy from
#' diseased individuals. Only applies for the method 'optimal.cutpoints.default').
#' @param tag.healthy	the value codifying healthy individuals in the status variable.
#' @param methods A character vector selecting the method(s) to be used: "CB" (cost-benefit method);
#' "MCT" (minimizes Misclassification Cost Term); "MinValueSp" (a minimum value set for Specificity);
#' "MinValueSe" (a minimum value set for Sensitivity); "ValueSe" (a value set for Sensitivity);
#' "MinValueSpSe" (a minimum value set for Specificity and Sensitivity); "MaxSp" (maximizes Specificity);
#' "MaxSe" (maximizes Sensitivity);"MaxSpSe" (maximizes Sensitivity and Specificity simultaneously);
#' "MaxProdSpSe" (maximizes the product of Sensitivity and Specificity or Accuracy Area);
#' "ROC01" (minimizes distance between ROC plot and point (0,1)); "SpEqualSe" (Sensitivity = Specificity);
#' "Youden" (Youden Index); "MaxEfficiency" (maximizes Efficiency or Accuracy, similar to minimize Error
#' Rate); "Minimax" (minimizes the most frequent error); "MaxDOR" (maximizes Diagnostic Odds Ratio);
#' "MaxKappa" (maximizes Kappa Index); "MinValueNPV" (a minimum value set for Negative Predictive Value);
#' "MinValuePPV" (a minimum value set for Positive Predictive Value); "ValueNPV" (a value set for Negative
#' Predictive Value);"ValuePPV" (a value set for Positive Predictive Value);"MinValueNPVPPV" (a minimum
#' value set for Predictive Values); "PROC01" (minimizes distance between PROC plot and point (0,1));
#' "NPVEqualPPV" (Negative Predictive Value = Positive Predictive Value); "MaxNPVPPV" (maximizes Positive
#' Predictive Value and Negative Predictive Value simultaneously); "MaxSumNPVPPV" (maximizes the sum of
#' the Predictive Values); "MaxProdNPVPPV" (maximizes the product of Predictive Values); "ValueDLR.Negative"
#' (a value set for Negative Diagnostic Likelihood Ratio); "ValueDLR.Positive" (a value set for Positive
#' Diagnostic Likelihood Ratio); "MinPvalue" (minimizes p-value associated with the statistical Chi-squared
#' test which measures the association between the marker and the binary result obtained on using the
#' cutpoint); "ObservedPrev" (The closest value to observed prevalence); "MeanPrev" (The closest value
#' to the mean of the diagnostic test values); or "PrevalenceMatching" (The value for which predicted
#' prevalence is practically equal to observed prevalence).
#' @param data A data frame containing all needed variables.
#' @param direction A Character string specifying the direction to compute the ROC curve. By default
#' individuals with a test value lower than the cutoff are classified as healthy (negative test),
#' whereas patients with a test value greater than (or equal to) the cutoff are classified as
#' diseased (positive test). If this is not the case, however, and the high values are related to
#' health, this argument should be established at ">".
#' @param categorical.cov A character string with the name of the categorical covariate according
#' to which optimal cutpoints are to be calculated. The default is NULL (no categorical covariate).
#' @param pop.prev The value of the disease's prevalence. The default is \code{NULL} (prevalence is
#' estimated on the basis of sample prevalence). It can be a vector indicating the prevalence values
#' for each categorical covariate level.
#' @param control Output of the \code{\link[OptimalCutpoints]{control.cutpoints}}  function.
#' @param ci.fit A logical value. If \code{TRUE}, inference is performed on the accuracy measures
#' at the optimal cutpoint. The default is \code{FALSE}.
#' @param conf.level A numerical value with the confidence level for the construction of the
#' confidence intervals. The default value is 0.95.
#' @param trace A logical value. If \code{TRUE}, information on progress is shown. The default is \code{FALSE}.
#' @param \code{...} Further arguments passed to or from other methods. None are used in this method.
#'
#' @details
#' See \code{\link[OptimalCutpoints]{optimal.cutpoints}}.
#'
#' @return
#' Returns an object of class "optimal.cutpoints" with the following components: \code{methods},
#' \code{levels.cat}, \code{call}, \code{dat}, etc. See \code{\link[OptimalCutpoints]{optimal.cutpoints}}
#' for details.
#'
#' @examples
#'    beta0 <- c(-2, 2, -2, 2, -2, -2)
#'    dat <- rdat(n=500, beta=beta0)
#'    fit <- glm(y ~., data=dat, family=binomial(link="logit"),
#'               control=glm.control(epsilon=1e-8, maxit=100, trace=FALSE))
#'    summary(fit)
#'    dat.test <- rdat(n=500, beta=-beta0)
#'    pred <- predict(fit, newdata=dat.test, type="response")
#'    dat.tmp <- data.frame(pred=pred, y=dat.test$y)
#'    result <- OptimalCut(X="pred", status="y", data=dat.tmp, tag.healthy=0, methods = "Youden")
#'    c0 <- result$output$Global$optimal.cutoff$cutoff;
#'    c0  # THE BEST CUTOFF POINT
#' @seealso \code{\link[OptimalCutpoints]{optimal.cutpoints}}
#' @references
#'\itemize{
#' \item Lopez-Raton, M., Rodriguez-Alvarez, M. X., Cadarso-Suarez, C., and Gude-Sampedro, F. (2014).
#' OptimalCutpoints: An R package for selecting optimal cutpoints in diagnostic tests. \emph{Journal of
#' Statistical Software}, \bold{61}(8), 1. URL \url{http://www.jstatsoft.org/v61/i08/}.
#' \item Zhang, Z., Su, X., et al. (2019+). Bootstrap Confidence Intervals for Optimal Cutoff Point
#' in Logistic Regression. Submitted.
#' }
#' @import OptimalCutpoints


OptimalThreshold <-  function(X, status=NULL, tag.healthy, methods, data, direction = c("<", ">"),
                     categorical.cov = NULL, pop.prev = NULL, control = control.cutpoints(),
                     ci.fit = FALSE, conf.level = 0.95, trace = FALSE, ...)
{
  if (missing(methods) || is.null(methods)) {
    stop("'methods' argument required.", call. = FALSE)
  }
  if (any(!(methods %in% c("CB", "MCT", "MinValueSp",
                           "MinValueSe", "ValueSp", "ValueSe",
                           "MinValueSpSe", "MaxSp", "MaxSe", "MaxSpSe",
                           "MaxProdSpSe", "ROC01", "SpEqualSe",
                           "Youden", "MaxEfficiency", "Minimax",
                           "MaxDOR", "MaxKappa", "MinValueNPV",
                           "MinValuePPV", "ValueNPV", "ValuePPV",
                           "MinValueNPVPPV", "PROC01", "NPVEqualPPV",
                           "MaxNPVPPV", "MaxSumNPVPPV", "MaxProdNPVPPV",
                           "ValueDLR.Negative", "ValueDLR.Positive",
                           "MinPvalue", "ObservedPrev", "MeanPrev",
                           "PrevalenceMatching")))) {
    stop("You have entered an invalid method.", call. = FALSE)
  }
  if (missing(data) || is.null(data)) {
    stop("'data' argument required.", call. = FALSE)
  }
  if (missing(X) || is.null(X)) {
    stop("'X' argument required.", call. = FALSE)
  }
  if (missing(status) || is.null(status)) {
    stop("'status' argument required.", call. = FALSE)
  }
  if (missing(tag.healthy) || is.null(tag.healthy)) {
    stop("'tag.healthy' argument required.", call. = FALSE)
  }
  if (is.logical(ci.fit) == FALSE) {
    stop("'ci.fit' must be a logical-type argument.",
         call. = FALSE)
  }
  if (conf.level < 0 | conf.level > 1 | length(conf.level) !=
      1) {
    stop("'conf.level' must be a single number between 0 and 1.",
         call. = FALSE)
  }
  if (is.logical(trace) == FALSE) {
    stop("'trace' must be a logical-type argument.",
         call. = FALSE)
  }
  if (is.null(pop.prev) & ci.fit == TRUE & !control$ci.PV %in%
      c("Exact", "Quadratic", "Wald", "AgrestiCoull",
        "RubinSchenker")) {
    warning(paste("Predictive Vaues CI: ``", control$ci.PV,
                  "'' method is not valid when prevalence is estimated from the sample.\n",
                  sep = ""), call. = FALSE)
  }
  if (!is.null(pop.prev) & ci.fit == TRUE & !control$ci.PV %in%
      c("Transformed", "NotTransformed", "GartNam")) {
    warning(paste("Predictive Values CI: \"", control$ci.PV,
                  "\" method is not valid when prevalence is not estimated from the sample.\n",
                  sep = ""), call. = FALSE)
  }
  direction <- match.arg(direction)
  if (!all(c(X, status, categorical.cov) %in% names(data))) {
    stop("Not all needed variables are supplied in 'data'.",
         call. = FALSE)
  }
  data <- na.omit(data[, c(X, status, categorical.cov)])
  res <- vector("list", length(methods))
  names(res) <- "output"     # PLACE REVISED
  if (!is.null(categorical.cov)) {
    if (!is.factor(data[, categorical.cov]))
      data[, categorical.cov] <- factor(data[, categorical.cov])
    data[, categorical.cov] <- droplevels(data[, categorical.cov])
    levels.cat <- levels(data[, categorical.cov])
    for (i in 1:length(methods)) {
      res[[i]] <- vector("list", length(levels.cat))
      names(res[[i]]) <- levels.cat
    }
  }
  else {
    levels.cat = 1
    res[[1]] <- vector("list", 1)
    names(res[[1]]) <- "Global"
  }
  pop.prev.new <- vector(length = length(levels(data[, categorical.cov])))
  if (is.null(pop.prev))
    pop.prev <- NA
  if (!is.null(categorical.cov) & length(pop.prev) != 1 & length(pop.prev) !=
      length(levels(data[, categorical.cov]))) {
    stop("You have entered different values for prevalence which \n do not coincide with categorical covariate levels.",
         call. = FALSE)
  }
  else if (!is.null(categorical.cov) & length(pop.prev) ==
           1) {
    pop.prev.new <- rep(pop.prev, length(levels(data[, categorical.cov])))
  }
  else if (is.null(categorical.cov) & length(pop.prev) > 1) {
    warning("You have entered several values for prevalence. \n The first value has been selected.",
            call. = FALSE, immediate. = TRUE)
    pop.prev.new <- pop.prev[1]
  }
  else {
    pop.prev.new <- pop.prev
  }
  for (i in 1:length(levels.cat)) {
    if (trace) {
      if (length(levels.cat) > 1) {
        text <- paste("Level: ", levels.cat[i],
                      sep = "")
        cat(text)
        cat("\nAnalysing ...\n\n")
      }
    }
    data.m <- if (length(levels.cat) != 1)
      data[data[, categorical.cov] == levels.cat[i], ]
    else data
    if (is.na(pop.prev.new[i])) {
      pop.prev.new[i] <- OptimalCutpoints:::calculate.sample.prev(data.m,
                                                                  status, tag.healthy)
    }
    OptimalCutpoints:::validate.prevalence(pop.prev.new[i])
    measures.acc <- OptimalCutpoints:::calculate.accuracy.measures(data = data.m,
                                                                   marker = X, status = status, tag.healthy = tag.healthy,
                                                                   direction = direction, pop.prev = pop.prev.new[i],
                                                                   control = control, conf.level = conf.level, ci.fit = ci.fit)
    for (j in 1:length(methods)) {
      if (trace) {
        text <- paste("Method: ", methods[j], sep = "")
        cat(text)
        cat("\nAnalysing ...\n\n")
      }
      res[[j]][[i]] <- eval(parse(text = paste("OptimalCutpoints:::function.",
                                               methods[j], sep = "")))(data = data.m,
                                                                       marker = X, status = status, tag.healthy = tag.healthy,
                                                                       direction = direction, pop.prev = pop.prev.new[i],
                                                                       control = control, conf.level = conf.level, ci.fit = ci.fit,
                                                                       measures.acc = measures.acc)
    }
  }
  res$methods <- methods
  if (length(levels.cat) != 1)
    res$levels.cat <- levels.cat
  res$call <- match.call()
  res$data <- data
  class(res) <- "optimal.cutpoints"
  invisible(res)
  res
}


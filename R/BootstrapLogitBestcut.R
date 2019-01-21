#' Bootstrapping the best cutoff points with logistic regression
#'
#' @name BootstrapLogitBestcut
#' @aliases BootstrapLogitOptimalCutpoint
#' @param formula An object of class \code{\link[stats]{formula}}, with the response on the left of
#' a \code{~} operator, and the terms on the right.
#' @param dat A data.frame in which to interpret the variables named in the \code{formula} argument.
#' @param B The number of bootstrap samples used. Suggest setting \eqn{B} large such as 1000 02 2000.
#' @param trace A logical value. If \code{TRUE}, information on progress is shown. The default is \code{FALSE}.
#' @param conf.level A numerical value with the confidence level for the construction of the
#' confidence intervals. The default value is 0.95.
#' @param cntr.glm The control setting for fitting logistic regression model with \code{glm}.
#' See \code{\link[stats]{glm.control}} for details.
#' @param seed The ranomd seed to reproduce results. Default is 123.
#' @param stratified.sampling A logical value indicating whether stratified sampling should be used in bootstrap
#' resampling and V-fold CV with stratification on the response. It is highly recommended with unbalanced
#' classification problems.
#' @param CV.method Values in \code{"none"}, \code{"OOB"}, \code{"LOO"}. Default is \code{"none"} with
#' no cross validation in computing the predicted probabilities. Otherwise, either \code{"OOB"} out-of-bag
#' samples or \code{"LOO"} leave-one-out jackknife method is used.
#' @param V.fold The number of folds in V-fold CV only when \code{CV.method="OOB"}. In this case, V-fold CV
#' is used to compute the predicted probabilities with the orinal sample data. The choice of \code{V.fold=3} is
#' suggested, as 1/3 = .33 roughly corresponds to .37 in OOB. Even so, the bootstrap/percentile CIs are still
#' doubtful.
#' @param cntr0 The control setting for \code{\link[OptimalCutpoints]{optimal.cutpoints}}, modified as
#' \code{OptimalThreshold}. See \code{\link[OptimalCutpoints]{control.cutpoints}} for details.
#' @param method.OptimalCutpoints A character vector selecting the method(s) to be used in \code{OptimalThreshold};
#' See \code{\link[OptimalCutpoints]{optimal.cutpoints}} for details. Choices are "CB" (cost-benefit method);
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
#'
#'#' @details
#' This immediate funciton runs a bootstrapping procedure for obtaining a bootstrap sample of the optimal cutpoint
#' with logistic regression. Its output will be used as input for \code{\link{BootstrapCI}}.
#'
#' @return
#' Four objects are retuned.
#' \describe{
#' \item{result.0}{The object from \code{OptimalThreshold}, same as \code{\link[OptimalCutpoints]{optimal.cutpoints}},
#' for the original sample data.}
#' \item{c.0}{The optimal cutoff point based on the original sample data;}
#' \item{c.bootstrap}{The bagging estimator for the optimal cutpoint;}
#' \item{N.Bn}{The incidence matrix of dim \eqn{latex}{B x n} associated with the bootstrap
#' resampling procedure, where element \eqn{latex}{B_{bi}} is the frequency that the \eqn{i}-th observation
#' shows up in the \eqn{b}-th bootstrap sample.}
#' }
#'
#' @examples
#'       data(fertility)
#'       form0 <- diagnosis ~  factor(season) + age + accident + fever + alcohol + sitting
#'       out.boots <- BootstrapLogitBestcut(formula=form0, dat=fertility, B=100,
#'                  trace=FALSE, conf.level=0.95, CV.method="none", V.fold=3,
#'                  cntr0=control.cutpoints(), method.OptimalCutpoints = "Youden")
#'       names(out.boots)
#'
#'@seealso \code{\link[OptimalCutpoints]{optimal.cutpoints}}, \code{\link{OptimalThreshold}}
#' @references
#'\itemize{
#' \item Zhang, Z., Su, X., et al. (2019+). Bootstrap Confidence Intervals for Optimal Cutoff Point
#' in Logistic Regression. Submitted.
#' }
#' @import OptimalCutpoints
#' @export


BootstrapLogitBestcut <- function(formula, dat, B=2000,
                                  trace=FALSE, conf.level=0.95,
                                  cntr.glm=glm.control(epsilon=1e-8, maxit=100, trace=FALSE),
                                  seed=123, stratified.sampling=FALSE,  # UNBALANCED CLASSIFICATION
                                  CV.method="none", V.fold=NULL,
                                  cntr0=control.cutpoints(), method.OptimalCutpoints = "Youden")
{
  set.seed(seed)
  form0 <- formula;
  n <- NROW(dat)
  yname <- all.vars(form0)[1]; ycol <- which(names(dat)==yname)
  y <- dat[, ycol]; y.levels <- sort(unique(y))
  if (length(y.levels)!=2) stop("Hmmm. A binary response is required!")
  if (min(mean(y==1), mean(y==0)) <= .1) print("For unbalanced classification, Set stratified.sampling=TRUE")
  if (stratified.sampling) {ID0 <- which(y==y.levels[1]); ID1 <- which(y==y.levels[2])}

  # COMPUTATION WITH THE ORIGINAL SAMPLE DATA
  # -------------------------------------------
  if (!is.element(CV.method, c("none", "LOO", "OOB"))) stop("Check CV.method, which must in none, OOB, or LOO!")
  if (CV.method=="LOO") {
    pred <- rep(0, n)
    for (i in 1:n) {
      fit.i <- glm(form0, data=dat[-i, ], family=binomial(link="logit"), control=cntr.glm)
      pred[i] <- predict(fit.i, newdata=dat[i,], type="response")
    }
  } else if (CV.method=="OOB" && !is.null(V.fold)) {
    print("In case of OOB, the bootstrap and percentile CI may not be valid.\n
          May just look at the IJ-based interval and its bias-corrected version.")
    if (!is.integer(V.fold) || V.fold <= 1) stop("V.fold needs to be an integer greater than 1. We suggest V.fold=3.")
    V <- V.fold
    if (stratified.sampling) {
      id.0 <- sample(1:V, size=length(ID0), replace=TRUE)
      id.1 <- sample(1:V, size=length(ID1), replace=TRUE)
      id.fold <- c(id.0, id.1)
    } else {id.fold <- sample(1:V, size=n, replace=TRUE)}
    subj <- pred <- NULL
    for (v in 1:V){
      train.v <- dat[id.fold!=v,];  test.v <- dat[id.fold==v,]
      fit.v <- glm(form0, data=train.v, family=binomial(link="logit"), control=cntr.glm)
      pred.v <- predict(fit.0, newdata=dat, type="response")
      subj.v <- which(id.fold==v)
      pred <- c(pred, pred.v); subj <- c(subj, subj.v)
    }
    pred <- pred[order(subj)]
  } else {
    fit.0 <- glm(form0, data=dat, family=binomial(link="logit"), control=cntr.glm)
    pred <- predict(fit.0, newdata=dat, type="response")
  }

  # OBTAIN THE OPTIMAL CUTPOINT ESTIMATE WITH SAMPLE DATA
  dat.tmp <- data.frame(pred=pred, y.obs=y)
  # result.0 <- optimal.cutpoints(pred~y.obs, data=dat.tmp, tag.healthy=0,
  # 	methods=method.OptimalCutpoints, control=cntr0)
  # c.0 <- result.0$Youden$Global$optimal.cutoff$cutoff
  # USING MODIFIED FUNCTION BestCut()
  # if (trace){print(table(dat.tmp$y.obs))}
  result.0 <- OptimalThreshold(X="pred", status="y.obs", data=dat.tmp, tag.healthy=0,
                               methods=method.OptimalCutpoints, control=cntr0, conf.level=conf.level)
  c.0 <- result.0$output$Global$optimal.cutoff$cutoff
  c.0 <- ifelse(length(c.0)>1, sample(c.0, 1), c.0)

  # START BOOTSTRAPPING
  # ---------------------
  c.bootstrap <- rep(0, B)
  N.Bn <- matrix(0, nrow=B, ncol=n)
  for (b in 1:B){
    if (trace){print(paste("-------------------- b = ", b, " ---------------------", sep=""))}
    # OBTAIN THE bTH BOOTSTRAP & OUT-OF-BAG SAMPLES
    if (stratified.sampling) {
      ID0 <- which(y==y.levels[1]); ID1 <- which(y==y.levels[2])
      id.0 <- sample(ID0, size=length(ID0), replace=TRUE)
      id.1 <- sample(ID1, size=length(ID1), replace=TRUE)
      id.b <- c(id.0, id.1)
    } else {id.b <- sample(x=1:n, size=n, replace=TRUE)}
    dat.b <- dat[id.b, ];
    if (CV.method=="OOB") {
      id.oob <- (1:n)[!is.element(1:n, sort(unique(id.b)))]
      dat.oob <- dat[id.oob, ]
    }
    N.Bn[b,] <- as.vector(table(factor(id.b, levels=1:n)))

    # FIT MODEL WITH D.b AND MAKE PREDICTION WITH D.oob
    fit.b <- glm(form0, data=dat.b, family=binomial(link="logit"),
                 control=cntr.glm)
    if (CV.method=="OOB") {
      pred <- predict(fit.b, newdata=dat.oob, type="response")
      y.obs <- dat.oob[, ycol]
    } else if (CV.method=="LOO") {
      pred <- rep(0, n)
      for (i in 1:n) {
        fit.i <- glm(form0, data=dat.b[-i, ], family=binomial(link="logit"), control=cntr.glm)
        pred[i] <- predict(fit.i, newdata=dat.b[i,], type="response")
      }
      y.obs <- dat.b[, ycol]
    } else {
      pred <- predict(fit.b, newdata=dat.b, type="response")
      y.obs <- dat.b[, ycol]
    }
    dat.tmp <- data.frame(pred=pred, y.obs=y.obs)
    if (trace){print(table(dat.tmp$pred, dat.tmp$y.obs))}
    # MODIFIED FUNCTION BestCut()
    result0 <- OptimalThreshold(X="pred", status="y.obs", data=dat.tmp, tag.healthy=0,
                                methods=method.OptimalCutpoints, control=cntr0, conf.level=conf.level)
    c0 <- result0$output$Global$optimal.cutoff$cutoff
    c.bootstrap[b] <- ifelse(length(c0)>1, sample(c0, 1), c0)
    if (trace) print(cbind(b=b, cutoff=c0[1]))
  }
  return(list(result.0=result.0, c.0=c.0, c.bootstrap=c.bootstrap, N.Bn=N.Bn))
}



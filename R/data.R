#' Fertility data concerning seminal quality
#' @docType data
#' @usage data(fertility)
#' @format The data are 100 observations on 10 variables, \code{diagnosis} being the target variable:
#' \describe{
#' \item{\code{season}}{Season in which the analysis was performed. 1) winter, 2) spring, 3) Summer, 4) fall.
#' Recoded as values (-1, -0.33, 0.33, 1) in the original data.}
#' \item{\code{age}}{Age at the time of analysis. 18-36. Standardized into (0, 1)}
#' \item{\code{childish}}{Childish diseases (ie , chicken pox, measles, mumps, polio): yes - 1 and no - 0.}
#' \item{\code{accident}}{Accident or serious trauma: yes (1) and no (0).}
#' \item{\code{surgical}}{Surgical intervention: Yes (1) and No (0).}
#' \item{\code{fever}}{High fevers in the last year: no (-1), less than three months ago (0), and more than three months ago (+1).}
#' \item{\code{alcohol}}{Frequency of alcohol consumption 1) several times a day, 2) every day,
#' 3) several times a week, 4) once a week, 5) hardly ever or never; Rescaled into continuous (0,1) values}
#' \item{\code{smoking}}{Smoking habit 1) never, 2) occasional 3) daily. Recoded as \{-1, 0, 1\}.}
#' \item{\code{sitting}}{Number of hours spent sitting per day. Rescaled into (0,1) values.}
#' \item{\code{diagnosis}}{Diagnosis normal (0), altered (1)}
#' }
#' @source
#' 100 volunteers provide a semen sample analyzed according to the WHO 2010 criteria. Sperm concentration
#' are related to socio-demographic data, environmental factors, health status, and life habits.
#'
#' @references
#'\itemize{
#' \item Gil, D., Girela, J. L., De Juan, J., Gomez-Torres, M. J., and Johnsson, M. (2012). Predicting seminal
#' quality with articial intelligence methods. \emph{Expert Systems with Applications}, \bold{39}(16): 12564--12573.
#' }
#' @examples
#'  data(fertility)
#'  head(fertility)
#'
"fertility"

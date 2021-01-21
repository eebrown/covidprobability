# False negative rate values reported by Kucirka et al.
# (point estimate and 95% CI)
# FNR_point <- c(100, 100, 97.7, 71, 38.7, 24.8, 20.1, 19.1, 20, 22.1, 25,
#                28.6, 32.5, 36.8, 41.2, 45.5, 49.6, 53.5, 57, 60.2, 63)/100
# FNR_lower <- c(100, 96, 58.8, 29.6, 18.4, 14, 12.5, 12, 12.8, 14.5, 16.7,
#                19.5, 22.6, 26.2, 30, 33.8, 37.6, 41.2, 44.6, 47.8, 50.8)/100
# FNR_upper <- c(100, 100, 99.9, 94.1, 64.6, 39.9, 31, 29.1, 30.2, 32.7, 36,
#                40.1, 44.4, 49, 53.6, 58, 62.1, 65.7, 68.9, 71.6, 74.2)/100
#
# sens_point <- (FNR_point - 1) * -1
# sens_upper <- (FNR_lower - 1) * -1
# sens_lower <- (FNR_upper - 1) * -1
#
# sens <- data.frame(point = sens_point, lower = sens_lower, upper = sens_upper)

#' COVID-19 PCR sensitivity by days since exposure
#'
#'
#'
#' @format A data frame with 21 rows and 3 variables:
#' \describe{
#'   \item{point}{point estimate of sensitivity}
#'   \item{lower}{lower 95% confidence interval of sensitivity}
#'   \item{upper}{upper 95% confidence interval of sensitivity}
#' }
#' @source \url{https://bmjopen.bmj.com/content/10/8/e039652}
"sens"


#' COVID-19 incubation period distribution
#'
#' \code{incubation <- data.frame(x = 0:21,y = dlnorm(c(0:21), 1.63, 0.5))}
#'
#' @format A data frame with 22 rows and 2 variables:
#' \describe{
#'   \item{x}{days since exposure}
#'   \item{y}{probability of symptoms emerging on this day}
#' }
#' @source \url{https://bmjopen.bmj.com/content/10/8/e039652}
"incubation"



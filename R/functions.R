
#' The remaining individuals who would not be expected to show symptoms yet
#'
#' Every day, a certain number of people are expected to show symptoms, based
#' on the incubation period. This would typically lead to further investigation
#' and ongoing suspicion of an outbreak. This function calculates the proportion
#' of individuals on a given day that would not be expected to have developed
#' symptoms yet. So if no one has developed symptoms, this proportion of people
#' could still have undetected COVID-19.
#'
#' @param t day
#' @param asympt The proportion of positive patients who would be expected not
#'               to ever develop symptoms (true asymptomatic patients).
#' @param mu The mean of a lognormal distribution that approximates the
#'           incubation period for COVID-19. E.g. 1.63 (see reference).
#' @param sigma The standard deviation of a lognormal distribution that
#'              approximates the incubation period for COVID-19. E.g. 0.5
#'              (see reference).
#'
#' @return Proportion who would not be expected to show symptoms yet
#' @references See McAloon et al. <https://bmjopen.bmj.com/content/10/8/e039652>
#' @export
#'
prop_remaining <- function(t, asympt, mu = 1.63, sigma = 0.5) {
  potential_sympt <- 1 - asympt
  prop_sympt_by <- stats::plnorm(t, mu, sigma)
  prop <- (potential_sympt - (potential_sympt * prop_sympt_by)) + asympt
  return(prop)
}


#' Calculate pretest probability change over time
#'
#' Calculates the pretest probability over time, assuming the individual does
#' not develop symptoms, by taking into account the distribution of incubation
#' periods (defined as the time from exposure to symptom onset).
#'
#' @param pre0 Initial pretest probability (on day of exposure)
#' @param asympt The proportion of positive patients who would be expected not
#'               to ever develop symptoms (true asymptomatic patients).
#' @param days The range of days (e.g. 0:14)
#' @param mu The mean of a lognormal distribution that approximates the
#'           incubation period for COVID-19. E.g. 1.63 (see reference).
#' @param sigma The standard deviation of a lognormal distribution that
#'              approximates the incubation period for COVID-19. E.g. 0.5
#'              (see reference).
#'
#' @return pretest probability by day (time series)
#' @references See McAloon et al. <https://bmjopen.bmj.com/content/10/8/e039652>
#' @export
#'
adjust_pretest <- function(pre0, asympt, days = 0:14,
                           mu = 1.63, sigma = 0.5) {
  return(pre0 * sapply(days, prop_remaining, asympt, mu, sigma))
}


#' Calculate posttest probability from pretest probability and test
#' characteristics
#'
#'
#'
#'
#' @param pretest_prob Pretest probability
#' @param sens Test sensitivity
#' @param spec Test specificity
#'
#' @return posttest probability
#' @export
#'
calc_postest_prob <- function(pretest_prob, sens, spec) {
  pretest_odds <- pretest_prob / (1 - pretest_prob)
  LR_neg <- (1 - sens) / spec
  posttest_odds <- pretest_odds * LR_neg
  postest_prob <- posttest_odds / (posttest_odds + 1)
  return(postest_prob)
}


#' Calculate post-test probability if testing occurred on each day in a series
#'
#' Given an initial pretest probability, and assuming symptoms never arise, with
#' each passing day the pretest probability will be lower, given the person did
#' not experience symptoms. This returns a vector of posttest probabilities
#' which takes all of the above into account, assuming a negative test on each
#' day. Note this is not a time series, and does not reflect if serial testing
#' were done each day and assumes testing was only done once.
#'
#' @param pre0 The pretest probability on day 0 (at exposure)
#' @param asympt The proportion of infected patients expected to remain
#'               asymptomatic throughout the course of infection
#' @param days The range of days (e.g. 0:14)
#' @param mu The mean of a lognormal distribution that approximates the
#'           incubation period for COVID-19. E.g. 1.63 (see reference).
#' @param sigma The standard deviation of a lognormal distribution that
#'              approximates the incubation period for COVID-19. E.g. 0.5
#'              (see reference).
#' @param sens A vector of sensitivities by day since exposure
#' @param spec The test specificity
#'
#' @return A vector of posttest probabilities
#' @export
#'
posttest_series <- function(pre0, asympt, days = 0:14, mu = 1.63, sigma = 0.5,
                            sens, spec) {
  pretest <- adjust_pretest(pre0 = pre0, asympt = asympt,
                            days = days, mu = mu, sigma = sigma)

  posttest <- data.frame(x=days, y=rep(NA, length(days)))

  for (t in seq_along(days)) {
    posttest$y[t] <- calc_postest_prob(pretest[t], sens[t], spec)
  }
  return(posttest)
}


#' Find the probability of any (at least one) event happening
#'
#' For an event that occurs with probability p, this function returns
#' the probability of an occurrence given n repetitions. p is numeric
#' and can be a vector.
#'
#' @param n The number of times to repeat the event (independent)
#' @param p The individual probability of the event happening
#' @return The probability of an event with the specified probability, after n
#'         repetitions
#' @export
#' @examples
#' probability_any(1, 0.5)
#' probability_any(2, 0.5)
#' probability_any(2, c(0.5, 1/3, 0.25))
probability_any <- function(n, p) {
  return(1 - ( (1 - p)^n))
}


#' Calculate a time series of probability for an individual following exposure
#'
#' The probability that an individual has COVID-19 will change over time as new
#' information is gleaned. The initial probability is the pretest probability
#' (pre0) associated with the nature of the interaction/exposure. This
#' probability will decrease with each passing day that the individual does not
#' develop symptoms. When a test is done, the probability is the posttest
#' probability; this reduces the probability based on the test characteristics
#' at the time of testing. Subsequently, the probability will continue to
#' decrease with each passing day that no symptoms develop. This function
#' returns a time series including those 3 phases.
#'
#' @param test_day Day of PCR test (days since exposure)
#' @param pre0 Pre-test probability of person on day of exposure
#' @param sens A vector of sensitivities by day since exposure
#' @param spec The specificity of the PCR test
#' @param asympt The proportion of infected patients expected to remain
#'               asymptomatic throughout the course of infection
#' @param days The range of days (e.g. 0:14)
#' @param mu The mean of a lognormal distribution that approximates the
#'           incubation period for COVID-19. E.g. 1.63 (see reference).
#' @param sigma The standard deviation of a lognormal distribution that
#'              approximates the incubation period for COVID-19. E.g. 0.5
#'              (see reference).
#' @return A time series of probabilities
#' @export
individual_probability <- function(test_day, pre0, sens, spec,
                                   asympt, days, mu, sigma) {

  pretest_series <- adjust_pretest(pre0, asympt,
                                   days, mu, sigma)[1:(test_day-1)]

  posttest <- posttest_series(pre0, asympt, days, mu, sigma, sens$point, spec)
  posttest_upper <- posttest_series(pre0, asympt, days,
                                    mu, sigma, sens$upper, spec)
  posttest_lower <- posttest_series(pre0, asympt, days,
                                    mu, sigma, sens$lower, spec)

  posttest_day <- posttest$y[test_day]
  posttest_day_upper <- posttest_upper$y[test_day]
  posttest_day_lower <- posttest_lower$y[test_day]

  after_series <- rep(posttest_day, utils::tail(days, 1)-test_day)
  after_series_upper <- rep(posttest_day_upper, utils::tail(days, 1)-test_day)
  after_series_lower <- rep(posttest_day_lower, utils::tail(days, 1)-test_day)

  additional_symptomatic <- c(sapply(days, prop_remaining, asympt, mu, sigma),
                              NA) - c(NA, sapply(days, prop_remaining, asympt,
                                                 mu, sigma))
  additional_symptomatic <-
    additional_symptomatic[(test_day+1):utils::tail(days, 1)] * -1

  after_series <- after_series -
                   (after_series * cumsum(additional_symptomatic))
  after_series_upper <- after_series_upper -
                   (after_series_upper * cumsum(additional_symptomatic))
  after_series_lower <- after_series_lower -
                   (after_series_lower * cumsum(additional_symptomatic))

  probability_series <- c(pretest_series, posttest_day, after_series)
  probability_series_upper <- c(rep(NA, test_day-1), posttest_day_upper,
                                after_series_upper)
  probability_series_lower <- c(rep(NA, test_day-1), posttest_day_lower,
                                after_series_lower)

  series <- data.frame(point = probability_series,
                       upper = probability_series_upper,
                       lower = probability_series_lower)
  return(series)
}

#' Calculate a time series of unit-wide probability following exposure
#'
#' To calculate the probability that any asymptomatic person has COVID-19,
#' this function treats each person/exposure as independent events and
#' calculates the probability time series using the individuals time series from
#' `individual_probability()`.
#'
#' @param test_day Day of PCR test (days since exposure)
#' @param pre0 Pre-test probability of person on day of exposure
#' @param sens A vector of sensitivities by day since exposure
#' @param spec The specificity of the PCR test
#' @param asympt The proportion of infected patients expected to remain
#'               asymptomatic throughout the course of infection
#' @param days The range of days (e.g. 0:14)
#' @param mu The mean of a lognormal distribution that approximates the
#'           incubation period for COVID-19. E.g. 1.63 (see reference).
#' @param sigma The standard deviation of a lognormal distribution that
#'              approximates the incubation period for COVID-19. E.g. 0.5
#'              (see reference).
#' @param n Number of exposed individuals
#' @return The probability of an event with the specified probability, after n
#'         repetitions
#' @export
unit_probability <- function(test_day, pre0, sens, spec, asympt, days, mu,
                             sigma, n) {

  series <- individual_probability(test_day = test_day, pre0 = pre0,
                                   sens = sens, spec = spec, asympt = asympt,
                                   days = days, mu = mu, sigma = sigma)

  series_n <- data.frame(point = probability_any(n, series$point),
                         upper = probability_any(n, series$upper),
                         lower = probability_any(n, series$lower))

  return(series_n)
}

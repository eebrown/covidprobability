
#' The remaining individuals who would not be expected to show symptoms yet
#'
#'
#'
#' @param t day
#' @param incubation The distribution of incubation period (days to symptom
#'                   onset) of the format
#'                   \code{data.frame(x = 0:21, y = dlnorm(c(0:21), 1.63, 0.5))}
#' @param asympt The proportion of positive patients who would be expected not
#'               to ever develop symptoms (true asymptomatic patients).
#'
#' @return Proportion who would not be expected to show symptoms yet
#' @export
#'
prop_remaining <- function(t, incubation, asympt) {
  potential_sympt <- 1 - asympt
  prop_sympt_by <- pracma::trapz(0:t, incubation$y[1:(t+1)])
  prop <- (potential_sympt - (potential_sympt * prop_sympt_by)) + asympt
  return(prop)
}



#' Calculate pretest probability change over time
#'
#' Calculates the pretest probability over time, assuming the individual does
#' not develop symptoms, by taking into account the distribution of incubation
#' periods (defined as the time from exposure to symptom onset).
#'
#' @param pretest Pretest probability (time series)
#' @param incubation The distribution of incubation period (days to symptom
#'                   onset) of the format
#'                   \code{data.frame(x = 0:21, y = dlnorm(c(0:21), 1.63, 0.5))}
#' @param asympt The proportion of positive patients who would be expected not
#'               to ever develop symptoms (true asymptomatic patients).
#'
#' @return posttest probability
#' @export
#'
adjust_pretest <- function(pretest, incubation, asympt) {
  return(pretest * sapply(incubation$x, prop_remaining, incubation, asympt))
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
#' @param incubation The distribution of incubation period (days to symptom
#'                   onset) of the format
#'                   \code{data.frame(x = 0:21, y = dlnorm(c(0:21), 1.63, 0.5))}
#' @param asympt The proportion of infected patients expected to remain
#'               asymptomatic throughout the course of infection
#' @param sens A vector of sensitivities by day since exposure
#' @param spec The test specificity
#'
#' @return A vector of posttest probabilities
#' @export
#'
posttest_series <- function(pre0, incubation, asympt, sens, spec) {
  pretest <- adjust_pretest(pre0, incubation, asympt)
  posttest <- data.frame(x=incubation$x, y=rep(NA, length(incubation$x)))
  for (t in 0:(length(incubation$x))) {
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
#'
#'
#' @param test_day Day of PCR test (days since exposure)
#' @param pre0 Pre-test probability of person on day of exposure
#' @param sens A vector of sensitivities by day since exposure
#' @param spec The specificity of the PCR test
#' @param incubation The distribution of incubation period (days to symptom
#'                   onset) of the format
#'                   \code{data.frame(x = 0:21, y = dlnorm(c(0:21), 1.63, 0.5))}
#' @param asympt The proportion of infected patients expected to remain
#'               asymptomatic throughout the course of infection
#' @return The probability of an event with the specified probability, after n
#'         repetitions
#' @export
individual_probability <- function(test_day, pre0, sens, spec, incubation,
                                   asympt) {

  pretest_series <- adjust_pretest(pre0, incubation, asympt)[1:(test_day-1)]

  posttest <- posttest_series(pre0, incubation, asympt, sens$point, spec)
  posttest_upper <- posttest_series(pre0, incubation, asympt, sens$upper, spec)
  posttest_lower <- posttest_series(pre0, incubation, asympt, sens$lower, spec)

  posttest_day <- posttest$y[test_day]
  posttest_day_upper <- posttest_upper$y[test_day]
  posttest_day_lower <- posttest_lower$y[test_day]

  after_series <- rep(posttest_day, 14-test_day)
  after_series_upper <- rep(posttest_day_upper, 14-test_day)
  after_series_lower <- rep(posttest_day_lower, 14-test_day)

  # additional proportion expected to have developed symptoms by this day
  prop_sympt_by <- function(t, dist=incubation) {
    pracma::trapz(0:t, dist$y[1:(t+1)])
  }

  additional_symptomatic <- c(sapply(0:14, prop_sympt_by), NA) -
                            c(NA, sapply(0:14, prop_sympt_by))
  additional_symptomatic <- additional_symptomatic[(test_day+1):14]

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
#'
#'
#' @param test_day Day of PCR test (days since exposure)
#' @param pre0 Pre-test probability of person on day of exposure
#' @param sens A vector of sensitivities by day since exposure
#' @param spec The specificity of the PCR test
#' @param incubation The distribution of incubation period (days to symptom
#'                   onset) of the format
#'                   \code{data.frame(x = 0:21, y = dlnorm(c(0:21), 1.63, 0.5))}
#' @param asympt The proportion of infected patients expected to remain
#'               asymptomatic throughout the course of infection
#' @param n Number of exposed individuals
#' @return The probability of an event with the specified probability, after n
#'         repetitions
#' @export
unit_probability <- function(test_day, pre0, sens, spec, incubation, asympt,
                             n) {

  series <- individual_probability(test_day = test_day, pre0 = pre0,
                                   sens = sens, spec = spec,
                                   incubation = incubation, asympt = asympt)

  series_n <- data.frame(point = probability_any(n, series$point),
                         upper = probability_any(n, series$upper),
                         lower = probability_any(n, series$lower))

  return(series_n)
}

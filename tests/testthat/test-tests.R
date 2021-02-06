test_that("prop_remaining produces expected values", {

  # On day 0, all should be remaining

  expect_equal(prop_remaining(t = 0, asympt = .28, mu = 1.63, sigma = 0.5),
               1)

  # After many days, should equal asymptomatic

  expect_equal(prop_remaining(t = 1000, asympt = .28, mu = 1.63, sigma = 0.5),
               0.28)


  # After the median incubation period, half of symptomatic cases done

  expect_equal(prop_remaining(t = exp(1.63), asympt = .28, mu = 1.63,
                              sigma = 0.5),
               0.28 + ((1-0.28)/2))

})

test_that("ajust_pretest runs", {

  # first set of values
  pre0 <- 0.6; asympt <- 0.2; days <- 21; mu <- 1.63; sigma <- 0.5

  a1 <- adjust_pretest(pre0 = pre0, asympt = asympt, days = days, mu = mu,
                       sigma = sigma)

  expect_true(all(a1 > asympt * pre0))
  expect_true(a1[1] < pre0)
  expect_equal(a1[5], pre0 * prop_remaining(5, asympt = asympt,
                                            mu = mu, sigma = sigma))

  #second set of values
  pre0 <- 0.2; asympt <- 0.05; days <- 20; mu <- 1.9; sigma <- 0.6

  a1 <- adjust_pretest(pre0 = pre0, asympt = asympt, days = days, mu = mu,
                       sigma = sigma)

  expect_true(all(a1 > asympt * pre0))
  expect_true(a1[1] < pre0)
  expect_equal(a1[5], pre0 * prop_remaining(5, asympt = asympt,
                                            mu = mu, sigma = sigma))

  #corner case
  pre0 <- 0.2; asympt <- 1; days <- 20; mu <- 1.9; sigma <- 0.6

  a1 <- adjust_pretest(pre0 = pre0, asympt = asympt, days = days, mu = mu,
                       sigma = sigma)

  expect_true(all(a1 == pre0))

  #corner case
  pre0 <- 1; asympt <- 0; days <- 200; mu <- 1.9; sigma <- 0.6

  a1 <- adjust_pretest(pre0 = pre0, asympt = asympt, days = days, mu = mu,
                       sigma = sigma)

  expect_true(a1[200] < 0.0001)

})


test_that("calc_postest_prob produces accurate values", {

  expect_equal(calc_postest_prob(pretest_prob = 0.9, sens = 0.8, spec = 0.8),
               0.6923077)

  expect_equal(calc_postest_prob(pretest_prob = 0.2, sens = 0.75, spec = 1),
               0.05882353)

})

test_that("posttest_series produces expected values", {

pre0 <- 0.2; asympt <- 0.28; days <- 14; mu <- 1.63; sigma <- 0.5; spec <- 1

test_posttest <- posttest_series(pre0 = pre0, asympt = asympt, days = days,
                                 mu = mu, sigma = sigma, sens = sens$point,
                                 spec = 1)

expect_equal(test_posttest$x[1], 1)

expect_equal(test_posttest$y[6],
             calc_postest_prob(adjust_pretest(pre0, asympt, days, mu, sigma)[6],
                               sens$point[6], 1))

})

test_that("probability_any produces accurate values", {

expect_equal(probability_any(4, 0.5), 1 - 1/16)

expect_equal(probability_any(2, 1/3), 5/9)

expect_equal(probability_any(4, 1/4), 1/4+3/16+9/64+27/256)

})

test_that("individual_probability produces expected values", {

test_day <- 9; pre0 <- 0.2; asympt <- 0.28; days <- 14; mu <- 1.63
sigma <- 0.5; spec <- 1

i1 <- individual_probability(test_day, pre0, sens, spec, asympt, days, mu,
                             sigma)

# expected value on the day of the test
test_day_post <- calc_postest_prob(pre0 * prop_remaining(t = test_day, asympt,
                                                         mu, sigma),
                              sens=sens$point[test_day], spec)

expect_equal(test_day_post, i1$point[test_day])

# expected value the day before the test
before_test_day <- adjust_pretest(pre0 = pre0, asympt = asympt, days = days,
                           mu = mu, sigma = sigma)

expect_equal(before_test_day[test_day-1], i1$point[test_day-1])
expect_equal(before_test_day[1:test_day-1], i1$point[1:test_day-1])


test_day <- 12; pre0 <- 0.3; asympt <- 0.5; days <- 17; mu <- 1.7; sigma <- 0.8
spec <- 0.9

i2 <- individual_probability(test_day, pre0, sens, spec, asympt, days, mu,
                             sigma)

test_day_post2 <- calc_postest_prob(pre0 * prop_remaining(t = test_day, asympt,
                                                          mu, sigma),
                                   sens=sens$point[test_day], spec)

expect_equal(test_day_post2, i2$point[test_day])

# expected value the day before the test
before_test_day2 <- adjust_pretest(pre0 = pre0, asympt = asympt, days = days,
                                  mu = mu, sigma = sigma)

expect_equal(before_test_day2[test_day-1], i2$point[test_day-1])
expect_equal(before_test_day2[1:test_day-1], i2$point[1:test_day-1])

# error with out of range days
expect_error(individual_probability(30, pre0, sens, spec, asympt, days, mu,
                                    sigma))

expect_error(individual_probability(7, pre0, sens, spec, asympt, 22, mu,
                                    sigma))


})

test_that("unit_probability produces expected values", {

  test_day <- 9; pre0 <- 0.2; asympt <- 0.28; days <- 14; mu <- 1.63
  sigma <- 0.5; spec <- 1; n <- 2

  i1 <- individual_probability(test_day, pre0, sens, spec, asympt, days, mu,
                               sigma)

  u1 <- unit_probability(test_day, pre0, sens, spec, asympt, days, mu,
                               sigma, n = n)

  expect_equal(probability_any(n, i1$point[4]), u1$point[4])

  n <- 30

  u1 <- unit_probability(test_day, pre0, sens, spec, asympt, days, mu,
                         sigma, n = n)

  expect_equal(probability_any(n, i1$point[9]), u1$point[9])

})

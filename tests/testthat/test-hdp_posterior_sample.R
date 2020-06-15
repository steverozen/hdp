
test_that("hdp_posterior_sample", {

  reg <- new.env()
  load("hdp.burnin.expected.Rdata", envir = reg)

  reg2 <- new.env()
  load("hdp.posterior.sample.expected.Rdata", envir = reg2)

  set.seed(44)

  retvalx <- hdp_posterior_sample(burnin.output = reg$retvalx,
                                  n = 15,
                                  space = 10,
                                  cpiter = 3,
                                  seed=44,
                                  verbosity=0)

  # save(retvalx, file = "hdp.posterior.sample.expected.Rdata")

  expect_equal(retvalx, reg2$retvalx)
})
